/// \file ijkcoord.txx
/// ijk templates for coordinate arithmetic
/// Version 0.4.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009-2017 Rephael Wenger

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

#ifndef _IJKCOORD_
#define _IJKCOORD_

#include <algorithm>
#include <cmath>

#include "ijk.txx"


namespace IJK {

  // **************************************************
  //! @name Basic coordinate operations
  // **************************************************

  ///@{

  /// Set all vertex coordinates to \a c.
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param c Scalar constant.  Set vertex coordinates to \a c.
  /// @param[out] coord[] Output coordinates.
  template <typename DTYPE, typename STYPE, typename CTYPE>
  void set_coord(const DTYPE dimension, const STYPE c, CTYPE coord[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord[d] = CTYPE(c); };
  }

  /// Copy \a coord0[] to \a coord1[].
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Input coordinates.
  /// @param[out] coord1[] Output coordinates.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void copy_coord(const DTYPE dimension, 
                  const CTYPE0 coord0[], CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = coord0[d]; };
  }

  /// Copy \a coord0[] to \a coord1[]. \n
  /// - Faster algorithm when \a coord0[] and \a coord1[] are both type float.
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Input coordinates.
  /// @param[out] coord1[] Output coordinates.
  template <typename DTYPE>
  void copy_coord(const DTYPE dimension, 
                  const float * coord0[], const float * coord1[])
  {
    std::copy(coord0, coord0+dimension, coord1);
  }

  /// Copy \a coord0[] to \a coord1[]. \n
  /// - Faster algorithm when \a coord0[] and \a coord1[] are both type double.
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Input coordinates.
  /// @param[out] coord1[] Output coordinates.
  template <typename DTYPE>
  void copy_coord(const DTYPE dimension, 
                  const double * coord0[], const double * coord1[])
  {
    std::copy(coord0, coord0+dimension, coord1);
  }

  /// Add \a coord0[] to \a coord1[].
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Input coordinates.
  /// @param coord1[] Input coordinates.
  /// @param[out] coord2[] Output coordinate.
  ///   - coord2[] = (\a coord0[] + \a coord1[]).
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_coord(const DTYPE dimension, const CTYPE0 coord0[],
                 const CTYPE1 coord1[], CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { coord2[d] = coord0[d] + coord1[d]; };
  }

  /// Subtract \a coord1[] from \a coord0[].
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param[out] coord2 = Output coordinate equal to 
  ///    (\a coord0[] - \a coord1[]).
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void subtract_coord
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord2[d] = coord0[d] - coord1[d]; };
  }

  /// Multiply \a coord0[] by the scalar \a s.
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param s = Scalar.
  /// @param[out] coord1 = Output coordinate equal to (\a s * \a coord0[]).
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void multiply_coord
  (const DTYPE dimension, const STYPE s, const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = s*coord0[d]; };
  }

  /// Divide \a coord0[] by the scalar \a s.
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param s = Non-zero scalar.  
  /// @pre Scalar s is not zero.
  /// @param[out] coord1 = Output coordinate equal to (\a s * \a coord0[]).
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void divide_coord
  (const DTYPE dimension, const STYPE s, const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = coord0[d]/s; };
  }

  /*!
   * Scale each \a coord0[d] by \a scale[d].
   * @param dimension Coordinate dimension (= number of coordinates.)
   * @param coord0[] Input coordinates.
   * @param scale[] Scale factors.
   * @param[out] coord1[] Output coordinates.
   *   - coord1[d] =  (\a scale[d] * \a coord0[d]).
   */
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void scale_coord
  (const DTYPE dimension, const STYPE scale[], const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = scale[d]*coord0[d]; };
  }

  /*!
   * Reverse scaling of each \a coord0[d].
   * @param dimension Coordinate dimension (= number of coordinates.)
   * @param coord0[] Input coordinates.
   * @param scale[] Scale factors.
   * @param[out] coord1[] Output coordinates.
   *   - coord1[d] =  (\a coord0[d] / \a scale[d]).
   * @pre \a scale[d] is not zero for any \a d.
   */
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void reverse_coord_scaling
  (const DTYPE dimension, const STYPE scale[], const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = coord0[d]/scale[d]; };
  }

  /*!
   *  Add \a s * \a coord0[] to \a coord1[]).
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param s Scaling factor.
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param[out] coord2 Output coordinate.
   *    - coord2[] = \a s* \a coord[0] + \a coord1[].
   */
  template <typename DTYPE, typename STYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_scaled_coord
  (const DTYPE dimension, const STYPE s, const CTYPE0 coord0[],
   const CTYPE1 coord1[], CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { coord2[d] = s*coord0[d] + coord1[d]; };
  }

  /// Compute absolute value of the each coord.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void abs_coord
  (const DTYPE dimension, const CTYPE0 coord0[], CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      coord1[d] = coord0[d];
      if (coord1[d] < 0) { coord1[d] = -coord1[d]; }
    }
  }

  /// Compute the midpoint of two coordinates.
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param[out] midpoint_coord = Coordinates of midpoint.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void compute_midpoint
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   CTYPE2 midpoint_coord[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { midpoint_coord[d] = (coord0[d] + coord1[d])/2.0; };
  }

  /// Compute sum of squares of coordinates.
  /// - Equivalent to computing the square of the magnitude.
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param[out] sum = sum of squares
  template <typename DTYPE, typename CTYPE0, typename STYPE>
  void compute_sum_of_squares(const DTYPE dimension, const CTYPE0 coord0[],
                              STYPE & sum)
  {
    sum = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { sum = sum + coord0[d]*coord0[d]; }
  }

  /// Compute magnitude of coordinate vector.
  /// - Equivalent to computing the square root of the sum 
  ///   of squares of the coordinates.
  /// - Fast implementation, but returns 0 when the value of (magnitude^2) 
  ///   is less than the machine epsilon of CTYPE0 or CTYPE1.
  /// @param dimension = Coordinate dimension (= number of coordinates.)
  /// @param coord0 = Input coordinates.
  /// @param[out] magnitude = magnitude
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void compute_magnitude(const DTYPE dimension, const CTYPE0 coord0[],
                         CTYPE1 & magnitude)
  {
    compute_sum_of_squares(dimension, coord0, magnitude);
    magnitude = std::sqrt(magnitude);
  }

  ///@}


  // **************************************************
  /// @name Vector operations
  // **************************************************

  ///@{

  /// Compute inner product of \a coord0[] and \a coord1[].
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Input coordinates.
  /// @param coord1[] Input coordinates.
  /// @param[out] product Inner product.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename STYPE>
  void compute_inner_product(const DTYPE dimension, const CTYPE0 coord0[],
                             const CTYPE1 coord1[], STYPE & product)
  {
    product = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { product = product + coord0[d]*coord1[d]; }
  }

  /*!
   *  Project vector \a v onto direction \a dir[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param v[] Input vector coordinates.
   *  @param dir[] Direction.
   *  @param[out] v_proj[] Projection of \a v onto \a dir[].
   *     - v_proj[] = inner_product(v[],dir[]) * dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2>
  void project_vector(const DTYPE dimension, const VTYPE0 v[], 
                      const VTYPE1 dir[], VTYPE2 v_proj[])
  {
    VTYPE0 inner_product;
    compute_inner_product(dimension, v, dir, inner_product);
    for (DTYPE d = 0; d < dimension; d++) 
      { v_proj[d] = dir[d]*inner_product; }
  }


  /*!
   *  Normalize vector.\n
   *  If magnitude of \a coord0[] is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of \a coord0[].
   *  @param[out] flag_zero True if \a coord1[] is set to the zero vector.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename MTYPE1>
  void normalize_vector
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE0 max_small_magnitude, CTYPE1 coord1[],
   MTYPE1 & magnitude, bool & flag_zero)
  { 
    compute_magnitude(dimension, coord0, magnitude); 
    if (magnitude > max_small_magnitude) {
      multiply_coord(dimension, 1.0/magnitude, coord0, coord1);
      flag_zero = false;
    }
    else {
      set_coord(dimension, 0, coord1);
      flag_zero = true;
    }
  }


  /*!
   *  Normalize vector.\n
   *  If magnitude of \a coord0[] is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set 
   *     to the zero vector. \n
   *  Version which does not return \a flag_zero.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of \a coord0[].
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename MTYPE1>
  void normalize_vector
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE0 max_small_magnitude, CTYPE1 coord1[],
   MTYPE1 & magnitude)
  { 
    compute_magnitude(dimension, coord0, magnitude); 
    if (magnitude > max_small_magnitude) {
      multiply_coord(dimension, 1.0/magnitude, coord0, coord1);
    }
    else {
      set_coord(dimension, 0, coord1);
    }
  }

  /*!
   *  Normalize vector.\n
   *  -  If magnitude of \a coord0[] is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set 
   *     to the zero vector. \n
   *  -  Version which does not return \a magnitude or \a flag_zero.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE, 
            typename CTYPE1>
  void normalize_vector
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE max_small_magnitude, CTYPE1 coord1[])
  { 
    double magnitude;
    normalize_vector
      (dimension, coord0, max_small_magnitude, coord1, magnitude);
  }

  /*!
   *  Compute magnitude of coordinate vector.\n
   *  -  Robust, more numerically stable, but more expensive, computation. \n
   *  If magnitude of \a coord0[] is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *    to \a max_small_magnitude, set coord1[] to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of \a coord0[].
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename MTYPE1>
  void normalize_vector_robust
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE0 max_small_magnitude, CTYPE1 coord1[], MTYPE1 & magnitude)
  { 
    MTYPE1 max_abs_coord;  // Use MTYPE1 in case CTYPE0 is int.

    magnitude = 0;
    if (dimension <= 0) { return; }

    // compute maximum absolute value of any coordinate.
    max_abs_coord = coord0[0];
    if (max_abs_coord < 0) { max_abs_coord = -max_abs_coord; };
    for (DTYPE d = 1; d < dimension; d++) {
      MTYPE1 abs_coord = coord0[d];
      if (abs_coord < 0) { abs_coord = -abs_coord; };
      if (max_abs_coord < abs_coord)
        { max_abs_coord = abs_coord; }
    }

    if (max_abs_coord > 0) {
      for (DTYPE d = 0; d < dimension; d++) {
        MTYPE1 c = coord0[d]/max_abs_coord;
        magnitude += c*c;
      }
      magnitude = std::sqrt(magnitude);
      magnitude = max_abs_coord*magnitude;
    }

    if (magnitude > max_small_magnitude) {
      divide_coord(dimension, magnitude, coord0, coord1);
    }
    else {
      set_coord(dimension, 0, coord1);
    }
  }


  /*!
   *  Compute unit vector from coord0[] to coord1[].
   *  -  If distance from \a coord0[] to \a coord1[]is less than or equal 
   *     to \a max_small_magnitude, then \a coord2[] is set to the zero vector.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param max_small_magnitude If distance from \a coord0[] to \a coord1[]
   *     is less than or equal to \a max_small_magnitude, then \a coord2[] 
   *      is set to the zero vector.
   *  @param[out] coord2[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of (\a coord1[] - \a coord0[]).
   *  @param[out] flag_zero True if \a coord2[] is set to the zero vector.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename CTYPE2, typename MTYPE1>
  void compute_unit_vector
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const MTYPE0 max_small_magnitude, CTYPE2 coord2[],
   MTYPE1 & magnitude, bool & flag_zero)
  {
    subtract_coord(dimension, coord1, coord0, coord2);
    normalize_vector
      (dimension, coord2, max_small_magnitude, coord2, magnitude, flag_zero);
  }


  /*!
   *  Compute component of vector orthogonal to given direction
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param v[] Input vector coordinates.
   *  @param dir[] Direction.
   *  @param[out] v_orth[] Component of \a v orthogonal to \a dir[].
   *     - v_orth[] = (\a v - v_proj[]) where v_proj[] is the projection
   *                  of \a v onto \a dir[].
   * @pre dir[] is a unit vector.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2>
  void compute_orthogonal_vector(const DTYPE dimension, const VTYPE0 v[], 
                                 const VTYPE1 dir[], VTYPE2 v_orth[])
  {
    VTYPE0 inner_product;
    compute_inner_product(dimension, v, dir, inner_product);
    for (DTYPE d = 0; d < dimension; d++) 
      { v_orth[d] = v[d]-dir[d]*inner_product; }
  }

  /// Compute three orthogonal unit vectors.
  /// - Return false if v0, v1 or v2 are (near) zero 
  ///   or v0 and v1 are collinear or v0, v1 and v2 are coplanar.
  /// - Note: max_small_magnitude is at or near zero can result
  ///   in extremely incorrect results.
  /// @param[out] fail_type Type of failure.
  /// - If fail_type == 0, then v0 is zero vector.
  /// - If fail_type == 1, then v0 and v1 are collinear.
  /// - If fail_type == 2, then v0, v1 and v2 are coplanar.
  /// - If function returns true, then fail_type is undefined.
  template <typename DTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE,
            typename WTYPE0, typename WTYPE1, typename WTYPE2,
            typename FTYPE>
  bool compute_orthogonal_unit_vectors
  (const DTYPE dimension, 
   const VTYPE0 v0[], const VTYPE1 v1[], const VTYPE2 v2[],
   const MTYPE max_small_magnitude,
   WTYPE0 w0[], WTYPE1 w1[], WTYPE2 w2[],
   FTYPE & fail_type)
  {
    double magnitude;
    double tolerance_for_check(0.001);
    double inner_product;

    set_coord(dimension, 0, w0);
    set_coord(dimension, 0, w1);
    set_coord(dimension, 0, w2);

    fail_type = 0;

    // Compute unit vector w0 = v0/|v0|.
    normalize_vector(dimension, v0, max_small_magnitude, w0, magnitude);
    if (magnitude <= max_small_magnitude) { return(false); }

    fail_type = 1;

    // Compute vector w1 orthogonal to w0.
    compute_orthogonal_vector(dimension, v1, w0, w1);
    normalize_vector(dimension, w1, max_small_magnitude, w1, magnitude);
    if (magnitude <= max_small_magnitude) { return(false); }

    // Check w1 is orthogonal to w0.
    compute_inner_product(dimension, w0, w1, inner_product);
    if (inner_product > tolerance_for_check || 
        inner_product < -tolerance_for_check) {
      set_coord(dimension, 0, w1);
      return(false);
    }

    fail_type = 2;

    // Compute vector w2 orthogonal to w0 and w1.
    compute_orthogonal_vector(dimension, v2, w0, w2);
    compute_orthogonal_vector(dimension, w2, w1, w2);
    normalize_vector(dimension, w2, max_small_magnitude, w2, magnitude);
    if (magnitude <= max_small_magnitude) { return(false); }

    // Check w2 is orthogonal to w1 and w0.
    compute_inner_product(dimension, w0, w2, inner_product);
    if (inner_product > tolerance_for_check || 
        inner_product < -tolerance_for_check) {
      set_coord(dimension, 0, w2);
      return(false);
    }
    compute_inner_product(dimension, w1, w2, inner_product);
    if (inner_product > tolerance_for_check || 
        inner_product < -tolerance_for_check) {
      set_coord(dimension, 0, w2);
      return(false);
    }

    return(true);
  }

  /// Compute three orthogonal unit vectors.
  /// - Does not set a fail type.
  template <typename DTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE,
            typename WTYPE0, typename WTYPE1, typename WTYPE2>
  bool compute_orthogonal_unit_vectors
  (const DTYPE dimension, 
   const VTYPE0 v0[], const VTYPE1 v1[], const VTYPE2 v2[],
   const MTYPE max_small_magnitude,
   WTYPE0 w0[], WTYPE1 w1[], WTYPE2 w2[])
  {
    int fail_type;

    bool flag =
      compute_orthogonal_unit_vectors
      (dimension, v0, v1, v2, max_small_magnitude, w0, w1, w2, fail_type);

    return(flag);
  }

  ///@}


  // **************************************************
  /// @name Compute angles
  // **************************************************

  ///@{

  /*!
   *  Compute cosine of the angle between two vectors.
   *  - If either vector has magnitude <= max_small_magnitude, 
   *    return zero and set flag_zero to true.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_angle Cosine of angle between coord0[] and coord1[].
   *  @param[out] flag_zero True, if coord0[] or coord1[]
   *    have zero magnitude. False, otherwise.
   *    - Note: Computed magnitude may be zero even if vectors are not zero.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_angle
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const MTYPE max_small_magnitude, COS_TYPE & cos_angle, bool & flag_zero)
  {
    CTYPE0 mag0;
    CTYPE1 mag1;
    compute_magnitude(dimension, coord0, mag0);
    compute_magnitude(dimension, coord1, mag1);
    if (mag0 > max_small_magnitude && mag1 > max_small_magnitude) {
      compute_inner_product(dimension, coord0, coord1, cos_angle);
      flag_zero = false;
      cos_angle = (cos_angle/mag0)/mag1;

      // Bound cos_angle to [-1,1].
      // Numerical error could cause cos_angle to be outside bounds.
      if (cos_angle < -1) { cos_angle = -1; }
      if (cos_angle > 1) { cos_angle = 1; }
    }
    else {
      cos_angle = 0;
      flag_zero = true;
    }
  }

  /*!
   *  Compute cosine of the angle between two vectors.
   *  - Version without max_small_magnitude parameter.
   *  - If either vector has (near) zero magnitude, return zero
   *    and set flag_zero to true.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename COS_TYPE>
  void compute_cos_angle
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   COS_TYPE & cos_angle, bool & flag_zero)
  {
    compute_cos_angle(dimension, coord0, coord1, 0.0, cos_angle, flag_zero);
  }

  /*!
   *  Compute the cosine of angle (coord0, coord1, coord2).
   *  - If either coord1 or coord2 have (nearly) the same coordinates
   *    as coord0, return 0 and set flag_identical to true.
   * @param[out] cos_angle Cosine of angle(coord0, coord1, coord2).
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename MTYPE>
  void compute_cos_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   const MTYPE max_small_magnitude, CTYPE3 & cos_angle, bool & flag_identical)
  {
    IJK::ARRAY<CTYPE0> u1(dimension);
    IJK::ARRAY<CTYPE1> u2(dimension);

    IJK::subtract_coord(dimension, coord0, coord1, u1.Ptr());
    IJK::subtract_coord(dimension, coord2, coord1, u2.Ptr());
    IJK::compute_cos_angle(dimension, u1.PtrConst(), u2.PtrConst(),
                           max_small_magnitude, cos_angle, flag_identical);
  }


  /*!
   *  Compute the cosine of angle (coord0, coord1, coord2).
   *  - Version without max_small_magnitude parameter.
   *  - If either coord1 or coord2 have (nearly) the same coordinates
   *    as coord0, return zero and set flag_identical to true.
   * @param[out] cos_angle Cosine of angle(coord0, coord1, coord2).
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3>
  void compute_cos_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   CTYPE3 & cos_angle, bool & flag_identical)
  {
    IJK::compute_cos_triangle_angle
      (dimension, coord0, coord1, coord2, 0.0, cos_angle, flag_identical);
  }


  /*!
   *  Compute the cosine of angle (iv0, iv1, iv2).
   * - Vertex coordinates are coord[iv0*dimension], coord[iv1*dimension],
   *   and coord[iv2*dimension].
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_triangle_angle_coord_list
  (const DTYPE dimension, const CTYPE * vertex_coord, 
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2, 
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_angle, bool & flag_identical)
  {
    const CTYPE * coord0 = vertex_coord+dimension*iv0;
    const CTYPE * coord1 = vertex_coord+dimension*iv1;
    const CTYPE * coord2 = vertex_coord+dimension*iv2;

    compute_cos_triangle_angle
      (dimension, coord0, coord1, coord2, max_small_magnitude, 
       cos_angle, flag_identical);
  }


  /*!
   *  Compute the cosine of angle (iv0, iv1, iv2).
   *  - Version without max_small_magnitude parameter.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename COS_TYPE>
  void compute_cos_triangle_angle_coord_list
  (const DTYPE dimension, const CTYPE * vertex_coord, 
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2, 
   COS_TYPE & cos_angle, bool & flag_identical)
  {
    compute_cos_triangle_angle_coord_list
      (dimension, vertex_coord, iv0, iv1, iv2, 0.0, cos_angle, flag_identical);
  }


  /*!
   *  Compute the cosine of smallest angle of triangle (coord0, coord1, coord2).
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if two or three triangles edges have
   *    length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename MTYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   const MTYPE max_small_magnitude, CTYPE3 & cos_min_angle, bool & flag_zero)
  {
    IJK::ARRAY<CTYPE0> u01(dimension);
    IJK::ARRAY<CTYPE1> u02(dimension);
    IJK::ARRAY<CTYPE2> u12(dimension);
    CTYPE0 mag01, mag02, mag12;
    CTYPE3 cos_angleA, cos_angleB, cos_angleC;

    IJK::subtract_coord(dimension, coord1, coord0, u01.Ptr());
    IJK::subtract_coord(dimension, coord2, coord0, u02.Ptr());
    IJK::subtract_coord(dimension, coord2, coord1, u12.Ptr());

    compute_magnitude(dimension, u01.PtrConst(), mag01);
    compute_magnitude(dimension, u02.PtrConst(), mag02);
    compute_magnitude(dimension, u12.PtrConst(), mag12);
    
    if (mag01 > max_small_magnitude && mag02 > max_small_magnitude &&
        mag12 > max_small_magnitude) {

      flag_zero = false;

      compute_inner_product
        (dimension, u01.PtrConst(), u02.PtrConst(), cos_angleA);
      cos_angleA = (cos_angleA/mag01)/mag02;

      // Note: Angle b is formed by vectors -u01 and u12.
      compute_inner_product
        (dimension, u01.PtrConst(), u12.PtrConst(), cos_angleB);
      cos_angleB = -(cos_angleB/mag01)/mag12;

      compute_inner_product
        (dimension, u02.PtrConst(), u12.PtrConst(), cos_angleC);
      cos_angleC = (cos_angleC/mag02)/mag12;

      if (cos_angleA > cos_angleB) { cos_min_angle = cos_angleA; }
      else { cos_min_angle = cos_angleB; }
      if (cos_min_angle < cos_angleC) { cos_min_angle = cos_angleC; }

      // Bound cos_min_angle to [-1,1] in case
      //   numerical error causes cos_min_angle to be outside bounds.
      if (cos_min_angle < -1) { cos_min_angle = -1; }
      if (cos_min_angle > 1) { cos_min_angle = 1; }
    }
    else if (mag01 > max_small_magnitude && mag02 > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u01.PtrConst(), u02.PtrConst(), cos_angleA);
      cos_min_angle = (cos_angleA/mag01)/mag02;
    }
    else if (mag01 > max_small_magnitude && mag12 > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u01.PtrConst(), u12.PtrConst(), cos_angleB);
      cos_min_angle = -(cos_angleB/mag01)/mag12;
    }
    else if (mag02 > max_small_magnitude && mag12 > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u02.PtrConst(), u12.PtrConst(), cos_angleB);
      cos_min_angle = (cos_angleC/mag02)/mag12;
    }
    else {
      flag_zero = true;
    }
  }


  /// Compute cosine of min/max of polygon angles of polygon.
  /// - If two adjacent polygon vertices have (near) identical coordinates,
  ///   the angle at those vertices is skipped.
  /// @param[out] cos_min Cosine of minimum polygon angle.
  ///    - Note: Min angle has maximimum cosine.
  /// @param[out] cos_max Cosine of maximum polygon angle.
  ///    - Note: Max angle has minimum cosine.
  /// @param[out] num_angle Number of cos angles computed (not skipped).
  template <typename DTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1,
            typename CTYPE, typename COS_TYPE0, typename COS_TYPE1>
  void compute_cos_min_max_polygon_angles
  (const DTYPE dimension, const VTYPE poly_vert[], const NTYPE0 num_poly_vert,
   const CTYPE vertex_coord[],
   COS_TYPE0 & cos_min, COS_TYPE1 & cos_max, NTYPE1 & num_angle)
  {
    cos_min = -1;
    cos_max = 1;
    num_angle = 0;

    if (num_poly_vert < 3) {
      // Ignore degenerate polygons with 1 or 2 vertices.
      return;
    }

    for (NTYPE0 i0 = 0; i0 < num_poly_vert; i0++) {
      NTYPE0 i1 = (i0+1)%num_poly_vert;
      NTYPE0 i2 = (i0+2)%num_poly_vert;
      VTYPE iv0 = poly_vert[i0];
      VTYPE iv1 = poly_vert[i1];
      VTYPE iv2 = poly_vert[i2];
      COS_TYPE0 cos_angle;
      bool flag_duplicate_point;
      IJK::compute_cos_triangle_angle_coord_list
        (dimension, vertex_coord, iv0, iv1, iv2, 
         cos_angle, flag_duplicate_point);

      if (!flag_duplicate_point) {
        num_angle++;
        if (cos_angle > cos_min) { cos_min = cos_angle; }
        if (cos_angle < cos_max) { cos_max = cos_angle; }
      }
    }
  }

  /// Compute min/max polygon angles of polgyon.
  /// - If two adjacent polygon vertices have (near) identical coordinates,
  ///   the angle at those vertices is skipped.
  /// @param[out] num_angle Number of angles computed (not skipped).
  template <typename DTYPE, typename VTYPE, typename NTYPE0, typename NTYPE1,
            typename CTYPE, typename ATYPE0, typename ATYPE1>
  void compute_min_max_polygon_angles
  (const DTYPE dimension, const VTYPE poly_vert[], const NTYPE0 num_poly_vert,
   const CTYPE vertex_coord[],
   ATYPE0 & min_angle, ATYPE1 & max_angle, NTYPE1 & num_angle)
  {
    ATYPE0 cos_min, cos_max;

    compute_cos_min_max_polygon_angles
      (dimension, poly_vert, num_poly_vert, vertex_coord,
       cos_min, cos_max, num_angle);
    min_angle = std::acos(cos_min) * 180.0/M_PI;
    max_angle = std::acos(cos_max) * 180.0/M_PI;
  }

  ///@}


  // **************************************************
  /// @name Compute dihedral angles
  // **************************************************

  ///@{

  /// Compute cosine of the dihedral angle between triangle (iv0, iv1, iv2)
  ///   and (iv0, iv1, iv3).
  /// @param[out] flag_collinear True if (iv0, iv1, iv2) or (iv0, iv1, iv3)
  ///   are (near) collinear.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3,
            typename COS_TYPE>
  void compute_cos_dihedral_angle
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2, const VTYPE3 iv3,
   COS_TYPE & cos_angle, bool & flag_fail)
  {
    const CTYPE * vertex0_coord = vertex_coord+dimension*iv0;
    const CTYPE * vertex1_coord = vertex_coord+dimension*iv1;
    const CTYPE * vertex2_coord = vertex_coord+dimension*iv2;
    const CTYPE * vertex3_coord = vertex_coord+dimension*iv3;
    CTYPE magnitude;
    IJK::ARRAY<CTYPE> w01(dimension);
    IJK::ARRAY<CTYPE> w02(dimension);
    IJK::ARRAY<CTYPE> w03(dimension);
    IJK::ARRAY<CTYPE> u02(dimension);
    IJK::ARRAY<CTYPE> u03(dimension);

    cos_angle = 1.0;

    compute_unit_vector
      (dimension, vertex0_coord, vertex1_coord, 0.0, 
       w01.Ptr(), magnitude, flag_fail);
    if (flag_fail) { return; }
    subtract_coord(dimension, vertex2_coord, vertex0_coord, w02.Ptr());
    subtract_coord(dimension, vertex3_coord, vertex0_coord, w03.Ptr());

    compute_orthogonal_vector(dimension, w02.Ptr(), w01.Ptr(), u02.Ptr());
    compute_orthogonal_vector(dimension, w03.Ptr(), w01.Ptr(), u03.Ptr());

    normalize_vector
      (dimension, u02.Ptr(), 0.0, u02.Ptr(), magnitude, flag_fail);
    if (flag_fail) { return; }
    normalize_vector
      (dimension, u03.Ptr(), 0.0, u03.Ptr(), magnitude, flag_fail);
    if (flag_fail) { return; }

    compute_inner_product(dimension, u02.Ptr(), u03.Ptr(), cos_angle);

    // Bound cos_angle to [-1,1].
    // Numerical error could cause cos_angle to be outside bounds.
    if (cos_angle < -1) { cos_angle = -1; }
    if (cos_angle > 1) { cos_angle = 1; }
  }


  /// Compute cosine of min/max dihedral angles of a tetrahedron.
  template <typename DTYPE, typename VTYPE, typename CTYPE,
            typename COS_TYPE0, typename COS_TYPE1, typename NUM_TYPE>
  void compute_cos_min_max_tetrahedron_dihedral_angles
  (const DTYPE dimension, const VTYPE tetrahedra_vert[], 
   const CTYPE * vertex_coord, COS_TYPE0 & cos_min, COS_TYPE1 & cos_max, 
   NUM_TYPE & num_angle)
  {
    const NUM_TYPE NUM_VERT_PER_TETRAHEDRON(4);

    cos_min = -1;
    cos_max = 1;
    num_angle = 0;

    for (int i0 = 0; i0 < NUM_VERT_PER_TETRAHEDRON; i0++) {
      const int iv0 = tetrahedra_vert[i0];

      for (int i1 = i0+1; i1 < NUM_VERT_PER_TETRAHEDRON; i1++) {
        const int iv1 = tetrahedra_vert[i1];

        int i2 = 0;
        while (i2 == i0 || i2 == i1) { i2++; }

        int i3 = 0;
        while (i3 == i0 || i3 == i1 || i3 == i2) { i3++; }

        const int iv2 = tetrahedra_vert[i2];
        const int iv3 = tetrahedra_vert[i3];

        // Compute dihedral angle between triangle (iv0, iv1, iv2) and
        //   triangle (iv0, iv1, iv3).

        double cos_angle;
        bool flag_collinear;
        compute_cos_dihedral_angle
          (dimension, vertex_coord, iv0, iv1, iv2, iv3, 
           cos_angle, flag_collinear);

        if (!flag_collinear) {
          num_angle++;
          if (cos_angle > cos_min) { cos_min = cos_angle; }
          if (cos_angle < cos_max) { cos_max = cos_angle; }
        }
      }
    }

  }

  ///@}


  // **************************************************
  /// @name Distance operations
  // **************************************************

  ///@{

  /// Compute square of distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE>
  void compute_distance_squared
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance_squared)
  {
    distance_squared = 0.0;
    for (DTYPE d = 0; d < dimension; d++) {
      CTYPE0 diff = coord0[d] - coord1[d];
      distance_squared += diff*diff;
    }
  }

  /// Compute distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE>
  void compute_distance
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance)
  {
    compute_distance_squared(dimension, coord0, coord1, distance);
    distance = std::sqrt(distance);
  }

  /// Compute square of scaled distance between two points.
  /// Scale by coordinate-wise divide.
  /// @pre scale[d]>0 for all d in range [0,dimension-1].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename SCALE_TYPE, typename DIST_TYPE>
  void compute_scaled_distance_squared_divide
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const SCALE_TYPE scale, DIST_TYPE & scaled_distance_squared)
  {
    scaled_distance_squared = 0.0;
    for (DTYPE d = 0; d < dimension; d++) {
      CTYPE0 diff = (coord0[d] - coord1[d])/scale[d];
      scaled_distance_squared += diff*diff;
    }
  }

  /// Compute scaled distance between two points.
  /// - Scale by coordinate-wise divide.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename SCALE_TYPE, typename DIST_TYPE>
  void compute_scaled_distance_divide
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const SCALE_TYPE scale, DIST_TYPE & scaled_distance)
  {
    compute_scaled_distance_squared_divide
      (dimension, coord0, coord1, scale, scaled_distance);
    scaled_distance = std::sqrt(scaled_distance);
  }

  /// Compute Linf distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE, typename AXIS_TYPE>
  void compute_Linf_distance
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance, AXIS_TYPE & axis)
  {
    DIST_TYPE x;

    distance = 0;
    axis = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d])
        { x = coord1[d] - coord0[d]; }
      else
        { x = coord0[d] - coord1[d]; }

      if (x > distance) {
        distance = x; 
        axis = d;
      }
    }
  }

  /// Compute Linf distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE>
  void compute_Linf_distance
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance)
  {
    DIST_TYPE x;

    distance = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d])
        { x = coord1[d] - coord0[d]; }
      else
        { x = coord0[d] - coord1[d]; }

      if (x > distance) { distance = x; }
    }
  }

  /*! 
   *  Compute signed (L2) distance from coord0[] to hyperplane through coord1[].
   *  @tparam DIST_TYPE Distance type.  DIST_TYPE should be a signed type.
   *  @param dimension Coordinate and vector dimensions.
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Hyperplane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to hyperplane.
   *  @param[out] distance Signed distance to hyperplane.
   *    - distance = inner product of (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_signed_distance_to_hyperplane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    // Compute inner product of (coord0[]-coord1[]) and orth_dir[].
    distance = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { distance = distance + (coord0[d]-coord1[d])*orth_dir[d]; }
  }

  /*!
   *  Compute (unsigned L2) distance from coord0[] 
   *     to hyperplane through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Hyperplane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to hyperplane.
   *  @param[out] distance Distance to hyperplane.
   *     - distance = the absolute value of inner product of 
   *                 (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_hyperplane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    VTYPE signed_distance;

    compute_signed_distance_to_hyperplane
      (dimension, coord0, coord1, orth_dir, signed_distance);

    if (signed_distance >= 0) 
      { distance = signed_distance; }
    else
      { distance = -signed_distance; }
  }

  /*!
   *  Compute (L2) distance squared from coord0[] to line through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance_squared
   *    Distance squared from coord0[] to line through coord1[].
   *    - distance_squared = Magnitude squared of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_squared_to_line
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance_squared)
  {
    // Compute inner product of (coord0[]-coord1[]) and dir[].
    DIST_TYPE x = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { x = x + (coord0[d]-coord1[d])*dir[d]; }

    // Compute magnitude squared of ((coord0[]-coord1[])-x*dir[]).
    distance_squared = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      DIST_TYPE y = (coord0[d]-coord1[d]) - x*dir[d];
      distance_squared = distance_squared + (y*y);
    }
  }

  /*!
   *  Compute (unsigned, L2) distance from coord0[] to line through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance Distance from coord0[] to line through coord1[].
   *    - distance = Magnitude of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_line
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance)
  {
    compute_distance_squared_to_line
      (dimension, coord0, coord1, dir, distance);
    distance = std::sqrt(distance);
  }

  ///@}


  // **************************************************
  /// @name Intersection operations
  // **************************************************  

  ///@{

  /// Compute the intersection of a line and a hyperplane.
  /// - Line contains \a coord0[] and has direction \a line_dir[].
  /// - Hyperplane contains \a coord1[] 
  ///   and has orthogonal direction \a orth_dir[].
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Coordinates of point contained by line. 
  /// @param line_dir[] Line direction.
  /// @param coord1[] Coordinates of point contained by hyperplane.
  /// @param orth_dir[] Direction orthgonal to \a line_dir[].
  /// @param min_abs_cos_angle_LO Minimum value of the absolute value
  ///     of the cosine of the angle between \a line_dir[] and \a orth_dir[].
  /// @param[out] intersection_point[] 
  ///               Intersection point of line and hyperplane.
  /// @param[out] flag_succeeded True if computation succeeded.
  ///    Computations succeeds if:
  ///      abs(inner_product(line_dir[],orth_dir[])) > min_abs_cos_angle_LO.
  /// @pre \a line_dir[] is a unit vector.
  /// @pre \a orth_dir[] is a unit vector.
  /// @pre \a min_abs_cos_angle_LO >= 0.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE0, typename VTYPE1, typename MTYPE,
            typename CTYPE2>
  void intersect_line_hyperplane
  (const DTYPE dimension, const CTYPE0 coord0[], const VTYPE0 line_dir[],
   const CTYPE1 coord1[], const VTYPE1 orth_dir[],
   const MTYPE min_abs_cos_angle_LO,
   CTYPE2 intersection_point[], bool & flag_succeeded)
  {
    VTYPE0 cos_angle;
    compute_inner_product(dimension, line_dir, orth_dir, cos_angle);

    if (std::abs(cos_angle) > min_abs_cos_angle_LO) {
      VTYPE0 distance, t;
      flag_succeeded = true;

      compute_signed_distance_to_hyperplane
        (dimension, coord0, coord1, orth_dir, distance);
      t = -distance/cos_angle;
      add_scaled_coord(dimension, t, line_dir, coord0, intersection_point);
    }
    else {
      flag_succeeded = false;
      IJK::copy_coord(dimension, coord1, intersection_point);
    }
  }

  ///@}


  // **************************************************
  /// @name Average operations
  // **************************************************  

  ///@{

  /*!
   *  Compute the weighted average of three coordinates.
   *  - Compute w0*coord0[]+w1*coord1[]+w2*coord2[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param w0 Weight of coord0[].
   *  @param coord0[] Input coordinates.
   *  @param w1 Weight of coord1[].
   *  @param coord1[] Input coordinates.
   *  @param w2 Weight of coord2[].
   *  @param coord2[] Input coordinates.
   *  @param[out] coord3 Output coordinate.
   *    - coord2[] = \a w0* \a coord[0] + \a w1*coord1[] + \a w2*coord2[].
   */
  template <typename DTYPE, 
            typename WTYPE0, typename WTYPE1, typename WTYPE2,
            typename CTYPE0, typename CTYPE1, typename CTYPE2, typename CTYPE3>
  void compute_weighted_average_3coord
  (const DTYPE dimension, 
   const WTYPE0 w0, const CTYPE0 coord0[],
   const WTYPE1 w1, const CTYPE1 coord1[], 
   const WTYPE2 w2, const CTYPE2 coord2[],
   CTYPE3 coord3[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { coord3[d] = w0*coord0[d] + w1*coord1[d] + w2*coord2[d]; };
  }
  
  ///@}

  // **************************************************
  /// @name Comparison operators
  // **************************************************

  ///@{

  /// Return true if coordinates are the same.
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Input coordinates.
  /// @param coord1[] Input coordinates.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_coord_equal(const DTYPE dimension, 
                      const CTYPE0 coord0[], const CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] != coord1[d]) { return(false); }; }

    return(true);
  }

  /// Return true if some coordinate is the same.
  /// @param dimension Coordinate dimension (= number of coordinates.)
  /// @param coord0[] Input coordinates.
  /// @param coord1[] Input coordinates.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_some_coord_equal(const DTYPE dimension, 
                           const CTYPE0 coord0[], const CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] == coord1[d]) { return(true); }; }

    return(false);
  }

  /// Return true if all elements of coord0[] are less than or equal
  ///   to all elements of coord1[].
  // *** SHOULD USE coord0[] and coord1[] IN PARAMETER LIST? ***
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_coord_less_than_or_equal
  (const DTYPE dimension, const CTYPE0 coord0, const CTYPE1 coord1)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] > coord1[d]) { return(false); }; };

    return(true);
  }

  /// Return true if coord0[] is contained in the rectangular region
  ///   with opposing corners coord1[] and coord2[].
  template <typename DTYPE, typename CTYPE0, 
            typename CTYPE1, typename CTYPE2>
  bool is_coord_in_rectangular_region
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d]) {
        if (coord0[d] < coord2[d]) { return(false); }
      }
      else if (coord0[d] > coord1[d]) {
        if (coord0[d] > coord2[d]) { return(false); }
      }
    }

    return(true);
  }

  /// Return true if coord0[] is contained in the interior 
  ///   of the rectangular region with opposing corners coord1[] and coord2[].
  template <typename DTYPE, typename CTYPE0, 
            typename CTYPE1, typename CTYPE2>
  bool is_coord_in_interior_of_rectangular_region
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d]) {
        if (coord0[d] <= coord2[d]) { return(false); }
      }
      else if (coord0[d] > coord1[d]) {
        if (coord0[d] >= coord2[d]) { return(false); }
      }
      else {
        // coord0[d] == coord1[d]
        return(false); 
      }
    }

    return(true);
  }


  /// Return true if scalar s0 is contained in the range [s1,s2],
  ///   i.e. Return true if (s1 <= s0 <= s2) or if (s1 >= s0 >= s2).
  template <typename STYPE0, typename STYPE1, typename STYPE2>
  bool is_scalar_in_range
  (const STYPE0 s0, const STYPE1 s1, const STYPE2 s2)
  {
    if (s0 < s1) {
      if (s0 < s2) { return(false); }
    }
    else if (s0 > s1) {
      if (s0 > s2) { return(false); }
    }

    return(true);
  }

  /// Return true if coord0[d] is in range [coord1[d], coord2[d]].
  /// Note: coord1[d] is not necessarily less than coord2[d].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  bool is_dth_coord_in_range
  (const DTYPE d, CTYPE0 coord0[], const CTYPE1 coord1[], 
   const CTYPE2 coord2[])
  {
    return(is_scalar_in_range(coord0[d], coord1[d], coord2[d]));
  }

  ///@}


  // **************************************************
  /// @name Rounding operators
  // **************************************************

  ///@{

  /// Round x to nearest integer.  
  template <typename T>
  T round(const T x)
  {
    return x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);
  }
  

  /// Round single coordinate
  template <int n, typename T>
  T round_coord(const T x)
  {
    T y = round(n*x);
    return(y/n);
  }

  /// Round single coordinate
  template <typename ITYPE, typename T>
  T round_coord(const ITYPE n, const T x)
  {
    T y = round(n*x);
    return(y/n);
  }

  /// Round coordinates.
  template <int n, typename DTYPE, typename T0, typename T1>
  void round_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round_coord<n>(v0[d]); }
  }

  /// Round coordinates.
  template <typename ITYPE, typename DTYPE, typename T0, typename T1>
  void round_coord
  (const ITYPE n, const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round_coord(n, v0[d]); }
  }

  // Specialized rounding operators.

  /// Round single coordinate to nearest 1/16'th.
  template <typename T>
  T round16_coord(const T x) { return(round_coord<16>(x)); }

  /// Round single coordinate to nearest 1/32'th.
  template <typename T>
  T round32_coord(const T x) { return(round_coord<32>(x)); }

  /// Round single coordinate to nearest 1/64'th.
  template <typename T>
  T round64_coord(const T x) { return(round_coord<64>(x)); }

  /// Round single coordinate to nearest 1/128'th.
  template <typename T>
  T round128_coord(const T x) { return(round_coord<128>(x)); }

  /// Round single coordinate to nearest 1/256'th.
  template <typename T>
  T round256_coord(const T x) { return(round_coord<256>(x)); }

  /// Round single coordinate to nearest 1/512'th.
  template <typename T>
  T round512_coord(const T x) { return(round_coord<512>(x)); }

  /// Round single coordinate to nearest 1/1024'th.
  template <typename T>
  T round1024_coord(const T x) { return(round_coord<1024>(x)); }

  /// Round coordinates to nearest 1/16'th.
  template <typename DTYPE, typename T0, typename T1>
  void round16_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round16_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/32'th.
  template <typename DTYPE, typename T0, typename T1>
  void round32_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round32_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/64'th.
  template <typename DTYPE, typename T0, typename T1>
  void round64_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round64_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/128'th.
  template <typename DTYPE, typename T0, typename T1>
  void round128_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round128_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/256'th.
  template <typename DTYPE, typename T0, typename T1>
  void round256_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round256_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/512'th.
  template <typename DTYPE, typename T0, typename T1>
  void round512_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round512_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/1024'th.
  template <typename DTYPE, typename T0, typename T1>
  void round1024_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round1024_coord(v0[d]); }
  }

  ///@}  

  // **************************************************
  /// @name Clamp operator
  // **************************************************

  ///@{

  /*!
   * Clamp x to interval [0,1].
   * - If (x < 0), then x = 0.
   * - If (x > 1), then x = 1.
   */
  template <typename T>
  T clamp01_coord(const T x)
  {
    if (x < 0) { return(0); }
    else if (x > 1) { return(1); }

    return(x);
  }

  ///@}

  // **************************************************
  // @name OLD DEPRECATED Comparison operators
  // **************************************************

  ///@{

  /// DEPRECATED: Replaced by is_coord_equal()
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool compare_coord(const DTYPE dimension, 
                     const CTYPE0 coord0[], const CTYPE1 coord1[])
  {
    return(is_coord_equal(dimension, coord0, coord1));
  }

  /// DEPRECATED: Replaced by is_coord_less_than_or_equal()
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_less_than_or_equal_to
  (const DTYPE dimension, const CTYPE0 coord0, const CTYPE1 coord1)
  {
    return(is_coord_less_than_or_equal(dimension, coord0, coord1));
  }

  ///@}


  // **************************************************
  /// @name 2D versions
  // **************************************************

  ///@{

  /// Subtract \a coord1[] from \a coord0[].
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param[out] coord2 = Output coordinate equal to (\a coord0[] - \a coord1[]).
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void subtract_coord_2D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 coord2[])
  { 
    const int DIM2 = 2;
    subtract_coord(DIM2, coord0, coord1, coord2); 
  }

  ///@}


  // **************************************************
  /// @name 3D versions
  // **************************************************

  ///@{

  /// Set all vertex coordinates to \a c.
  /// @param c = Scalar constant.  Set vertex coordinates to \a c.
  /// @param[out] coord[] = Output coordinates.
  template <typename STYPE, typename CTYPE>
  void set_coord_3D(const STYPE c, CTYPE coord[])
  { 
    const int DIM3 = 3;
    set_coord(DIM3, c, coord); 
  }

  /// Copy \a coord0[] to \a coord1[].
  /// @param coord0 = Input coordinates.
  /// @param[out] coord1 = Output coordinates.
  template <typename CTYPE0, typename CTYPE1>
  void copy_coord_3D(const CTYPE0 coord0[], CTYPE1 coord1[])
  { 
    const int DIM3 = 3;
    copy_coord(DIM3, coord0, coord1); 
  }

  /// Add \a coord0[] to \a coord1[].
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param[out] coord2 = Output coordinate equal to (\a coord0[] + \a coord1[]).
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_coord_3D(const CTYPE0 coord0[],
                    const CTYPE1 coord1[], CTYPE2 coord2[])
  { 
    const int DIM3 = 3;
    add_coord(DIM3, coord0, coord1, coord2); 
  }

  /// Subtract \a coord1[] from \a coord0[].
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param[out] coord2 = Output coordinate equal to (\a coord0[] - \a coord1[]).
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void subtract_coord_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 coord2[])
  { 
    const int DIM3 = 3;
    subtract_coord(DIM3, coord0, coord1, coord2); 
  }

  /// Multiply \a coord0[] by the scalar \a s.
  template <typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void multiply_coord_3D
  (const STYPE s, const CTYPE0 coord0[], CTYPE1 coord1[])
  { 
    const int DIM3 = 3;
    multiply_coord(DIM3, s, coord0, coord1); 
  }

  /*!
   *  Add \a s * \a coord0[] to \a coord1[]).
   *  @param s Scaling factor.
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param[out] coord2[] Output coordinate.
   *     - coord2[] = (\a s * \a coord0[] + \a coord1[]).
   */
  template <typename STYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_scaled_coord_3D
  (const STYPE s, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 coord2[])
  {
    const int DIM3 = 3;
    add_scaled_coord(DIM3, s, coord0, coord1, coord2);
  }

  /// Compute absolute value of the each coord.
  template <typename CTYPE0, typename CTYPE1>
  void abs_coord_3D(const CTYPE0 coord0[], CTYPE1 coord1[])
  {
    const int DIM3 = 3;
    abs_coord(DIM3, coord0, coord1);
  }

  /// Compute the midpoint of two coordinates.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void compute_midpoint_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 midpoint_coord[])
  { 
    const int DIM3 = 3;
    compute_midpoint(DIM3, coord0, coord1, midpoint_coord); 
  }

  /// Compute inner product of \a coord0[] and \a coord1[].
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  /// @param[out] product = Inner product.
  template <typename CTYPE0, typename CTYPE1, typename STYPE>
  void compute_inner_product_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], STYPE & product)
  { 
    const int DIM3 = 3;
    compute_inner_product(DIM3, coord0, coord1, product); 
  }

  /// Compute sum of squares of coordinates.
  template <typename CTYPE0, typename STYPE>
  void compute_sum_of_squares_3D(const CTYPE0 coord0[], STYPE & sum)
  { 
    const int DIM3 = 3;
    compute_sum_of_squares(DIM3, coord0, sum); 
  }

  /// Compute magnitude of coordinate vector.
  /// @param coord0 = Input coordinates.
  /// @param[out] magnitude = magnitude
  template <typename CTYPE0, typename STYPE>
  void compute_magnitude_3D(const CTYPE0 coord0[], STYPE & magnitude)
  { 
    const int DIM3 = 3;
    compute_magnitude(DIM3, coord0, magnitude); 
  }

  /// Return true if coordinates are the same.
  /// @param coord0 = Input coordinates.
  /// @param coord1 = Input coordinates.
  template <typename CTYPE0, typename CTYPE1>
  bool is_coord_equal_3D(const CTYPE0 coord0[], const CTYPE1 coord1[])
  { 
    const int DIM3 = 3;
    return(is_coord_equal(DIM3, coord0, coord1)); 
  }

  /// Return true if all elements of coord0[] are less than or equal
  ///   to all elements of coord1[].
  template <typename CTYPE0, typename CTYPE1>
  bool is_coord_less_than_or_equal_3D
  (const CTYPE0 coord0, const CTYPE1 coord1)
  { 
    const int DIM3 = 3;
    return(is_coord_less_than_or_equal(DIM3, coord0, coord1)); 
  }

  //@}


  // **************************************************
  /// @name 2D and 3D Matrix operations
  // **************************************************

  ///@{

  /// Return determinant of a 2x2 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_2x2
  (const COORD_TYPE a00, const COORD_TYPE a01,
   const COORD_TYPE a10, const COORD_TYPE a11,
   RESULT_TYPE & result)
  {
    result = a00*a11 - a01*a10;
  }

  /// Return determinant of a 3x3 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_3x3
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], RESULT_TYPE & result)
  {
    RESULT_TYPE D;
    determinant_2x2(p1[1], p1[2], p2[1], p2[2], D);
    result = p0[0] * D;
    determinant_2x2(p2[1], p2[2], p0[1], p0[2], D);
    result += p1[0] * D;
    determinant_2x2(p0[1], p0[2], p1[1], p1[2], D);
    result += p2[0] * D;
  }


  /// Return the square of the Froebenius norm of a 3x3 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void compute_froebenius_norm_squared_3x3
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], RESULT_TYPE & result)
  {
    RESULT_TYPE x;
    compute_sum_of_squares_3D(p0, result);
    compute_sum_of_squares_3D(p1, x);
    result += x;
    compute_sum_of_squares_3D(p2, x);
    result += x;
  }


  /// Return the Froebenius norm of a 3x3 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void compute_froebenius_norm_3x3
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], RESULT_TYPE & result)
  {
    compute_froebenius_norm_squared_3x3(p0, p1, p2, result);
    result = std::sqrt(result);
  }


  /// Transpose a 3x3 matrix.
  template <typename COORD_TYPE>
  void transpose_3x3
  (COORD_TYPE p0[3], COORD_TYPE p1[2], COORD_TYPE p2[3])
  {
    std::swap(p0[1], p1[0]);
    std::swap(p0[2], p2[0]);
    std::swap(p1[2], p2[1]);
  }


  /// Return the first row of the co-factor matrix of a 3x3 matrix.
  template <typename COORD_TYPEA, typename COORD_TYPEB>
  void co_factor_3x3_row1
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3], COORD_TYPEB q0[3])
  {
    const int DIM3 = 3;

    for (int j0 = 0; j0 < DIM3; j0++) {
      const int j1 = (j0+1)%DIM3;
      const int j2 = (j0+2)%DIM3;
      determinant_2x2(p1[j1], p1[j2], p2[j1], p2[j2], q0[j0]);
    }
  }


  /// Return the co-factor matrix of a 3x3 matrix.
  template <typename COORD_TYPEA, typename COORD_TYPEB>
  void co_factor_3x3
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3],
   COORD_TYPEB q0[3], COORD_TYPEB q1[2], COORD_TYPEB q2[3])
  {
    co_factor_3x3_row1(p0, p1, p2, q0);
    co_factor_3x3_row1(p1, p2, p0, q1);
    co_factor_3x3_row1(p2, p0, p1, q2);
  }


  /// Return the transpose of the co-factor matrix of a 3x3 matrix.
  template <typename COORD_TYPEA, typename COORD_TYPEB>
  void transpose_co_factor_3x3
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3],
   COORD_TYPEB q0[3], COORD_TYPEB q1[2], COORD_TYPEB q2[3])
  {
    co_factor_3x3(p0, p1, p2, q0, q1, q2);
    transpose_3x3(q0, q1, q2);
  }


  /*!
   *  Compute inverse of a 3x3 matrix.
   *  If absolute values of the determinant is less than or equal to
   *     max_small_determinant, then result is set to the transpose
   *     of the co-factor matrix.
   *  @param[out] determinant Determinant.
   *  @param[out] flag_zero True if determinant is less
   *     than or equal to max_small_determinant.
   *  @pre max_small_determinant >= 0.
   */
  template <typename COORD_TYPEA, typename COORD_TYPEB,
            typename MTYPE, typename DET_TYPE>
  void compute_inverse_3x3
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3],
   const MTYPE max_small_determinant, 
   COORD_TYPEB q0[3], COORD_TYPEB q1[2], COORD_TYPEB q2[3],
   DET_TYPE & determinant, bool & flag_zero)
  { 
    const int DIM3(3);
    DET_TYPE D;

    // Initialize flag_zero.
    flag_zero = false;

    transpose_co_factor_3x3(p0, p1, p2, q0, q1, q2);
    determinant_3x3(p0, p1, p2, determinant);
    D = determinant;
    if (D < 0) { D = -D; }

    if (D <= max_small_determinant) {
      flag_zero = true;
    }
    else {
      divide_coord(DIM3, determinant, q0, q0);
      divide_coord(DIM3, determinant, q1, q1);
      divide_coord(DIM3, determinant, q2, q2);
    }
  }


  /// Return determinant of three 2D points. <br>
  /// More precisely, return the determinant 
  ///   of the 2x2 matrix [(p1-p0),(p2-p0)].
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_point_2D
  (const COORD_TYPE p0[2], const COORD_TYPE p1[2], const COORD_TYPE p2[2],
   RESULT_TYPE & result)
  {
    RESULT_TYPE v1[2], v2[2];

    subtract_coord_2D(p1, p0, v1);
    subtract_coord_2D(p2, p0, v2);

    determinant_2x2(v1[0], v1[1], v2[0], v2[1], result);
  }

  /// Return determinant of four 3D points. <br>
  /// More precisely, return the determinant 
  ///   of the 3x3 matrix [(p1-p0),(p2-p0),(p3-p0)]. <br>
  /// Note: Equivalent to -determinant((p0,1);(p1,1);(p2,1);(p3,1)).
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_point_3D
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], const COORD_TYPE p3[3],
   RESULT_TYPE & result)
  {
    RESULT_TYPE v1[3], v2[3], v3[3];

    subtract_coord_3D(p1, p0, v1);
    subtract_coord_3D(p2, p0, v2);
    subtract_coord_3D(p3, p0, v3);

    determinant_3x3(v1, v2, v3, result);
  }

  ///@}


  // **************************************************
  /// @name 3D vector operations
  // **************************************************

  ///@{

  /// Compute 3D unit vector.
  template <typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename CTYPE2, typename MTYPE1>
  void compute_unit_vector_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const MTYPE0 max_small_magnitude, CTYPE2 coord2[],
   MTYPE1 & magnitude, bool & flag_zero)
  {
    const int DIM3 = 3;
    compute_unit_vector(DIM3, coord0, coord1, max_small_magnitude, coord2,
                        magnitude, flag_zero);
  }

  ///@}


  // **************************************************
  /// @name 3D distances
  // **************************************************

  ///@{

  /*! 
   *  Compute signed (L2) distance from coord0[] to hyperplane through coord1[].
   *  @tparam DIST_TYPE Distance type.  DIST_TYPE should be a signed type.
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Plane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to hyperplane.
   *  @param[out] distance Signed distance to hyperplane.
   *    - distance = inner product of (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_signed_distance_to_plane_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    const int DIM3 = 3;
    compute_signed_distance_to_hyperplane
      (DIM3, coord0, coord1, orth_dir, distance);
  }

  /*!
   *  Compute (unsigned, L2) distance from coord0[] to plane through coord1[].
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Plane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to hyperplane.
   *  @param[out] distance Distance to hyperplane.
   *     - distance = the absolute value of inner product of 
   *                 (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_plane_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    const int DIM3 = 3;
    compute_distance_to_hyperplane
      (DIM3, coord0, coord1, orth_dir, distance);
  }

  /*!
   *  Compute (L2) distance squared from coord0[] to line through coord1[].
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance_squared
   *    Distance squared from coord0[] to line through coord1[].
   *    - distance_squared = Magnitude squared of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_squared_to_line_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance_squared)
  {
    const int DIM3 = 3;
    compute_distance_squared_to_line
      (DIM3, coord0, coord1, dir, distance_squared);
  }

  /*!
   *  Compute (unsigned, L2) distance from coord0[] to line through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance Distance from coord0[] to line through coord1[].
   *    - distance = Magnitude of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_line_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance)
  {
    const int DIM3 = 3;
    compute_distance_to_line(DIM3, coord0, coord1, dir, distance);
  }

  ///@}


  // **************************************************
  /// @name 3D intersections
  // **************************************************

  ///@{

  /// Compute the intersection of a line and a plane.
  /// Line contains \a coord0[] and has direction \a line_dir[].
  /// Plane contains \a coord1[] and has orthogonal direction \a orth_dir[].
  /// @param coord0[] Coordinates of point contained by line. 
  /// @param line_dir[] Line direction.
  /// @param coord1[] Coordinates of point contained by hyperplane.
  /// @param orth_dir[] Direction orthgonal to \a line_dir[].
  /// @param min_abs_cos_angle_LO Minimum value of the absolute value
  ///     of the cosine of the angle between \a line_dir[] and \a orth_dir[].
  /// @param[out] intersection_point[] 
  ///               Intersection point of line and hyperplane.
  /// @param[out] flag_succeeded True if computation succeeded.
  ///    Computations succeeds if:
  ///      abs(inner_product(line_dir[],orth_dir[])) > min_abs_cos_angle_LO.
  /// @pre \a line_dir[] is a unit vector.
  /// @pre \a orth_dir[] is a unit vector.
  /// @pre \a min_abs_cos_angle_LO >= 0.
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE0, typename VTYPE1, typename MTYPE,
            typename CTYPE2>
  void intersect_line_plane_3D
  (const CTYPE0 coord0[], const VTYPE0 line_dir[],
   const CTYPE1 coord1[], const VTYPE1 orth_dir[],
   const MTYPE min_abs_cos_angle_LO,
   CTYPE2 intersection_point[], bool & flag_succeeded)
  {
    const int DIM3 = 3;
    intersect_line_hyperplane
      (DIM3, coord0, line_dir, coord1, orth_dir, min_abs_cos_angle_LO,
       intersection_point, flag_succeeded);
  }

  ///@}


  // **************************************************
  /// @name 3D cross product
  // **************************************************

  ///@{

  /// Compute the cross product of two 3D vectors.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void compute_cross_product_3D
  (const VTYPE0 v0[], const VTYPE1 v1[], VTYPE2 v2[])
  {
    determinant_2x2(v0[1], v0[2], v1[1], v1[2], v2[0]);
    determinant_2x2(v0[0], v0[2], v1[0], v1[2], v2[1]);
    v2[1] = -v2[1];
    determinant_2x2(v0[0], v0[1], v1[0], v1[1], v2[2]);
  }

  ///@}

}

#endif
