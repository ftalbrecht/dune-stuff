// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Sven Kaulmann

#ifndef DUNE_STUFF_GRID_INTERSECTION_HH
#define DUNE_STUFF_GRID_INTERSECTION_HH

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#if HAVE_DUNE_GRID
# include <dune/grid/common/gridview.hh>
#endif

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/type_utils.hh>

namespace Dune {
namespace Stuff {
namespace Grid {


#if HAVE_DUNE_GRID


template< class GridPartOrViewType >
class Intersection
{
  template< class GridViewType, bool is_view >
  struct Choose
  {
    typedef typename GridViewType::Intersection Type;
  };

  template< class GridPartType >
  struct Choose< GridPartType, false >
  {
    typedef typename GridPartType::IntersectionType Type;
  };

  static const bool this_is_a_grid_view =
      std::is_base_of< GridView< typename GridPartOrViewType::Traits >, GridPartOrViewType >::value;

public:
  typedef typename Choose< GridPartOrViewType, this_is_a_grid_view >::Type Type;
}; // class Intersection


#endif // HAVE_DUNE_GRID


/**
  \brief      prints some basic information about a Dune::Intersection, namely the number of its corners and the
              coordinates of those corners.
  \tparam     IntersectionType
              Dune::Intersection compatible
  \param[in]  intersection
              Dune::Intersection, whose information should be printed
  \param[out] stream
              std::ostream, into which the information is printed
  **/
template< class IntersectionType >
void printIntersection(const IntersectionType& intersection,
                       std::ostream& out = std::cout,
                       const std::string prefix = "")
{
  out << prefix << Common::Typename< IntersectionType >::value() << std::endl;
  const auto& geometry = intersection.geometry();
  for (auto ii : DSC::valueRange(geometry.corners()))
    out << prefix
        << "  corner " + Common::toString(ii) << " = " << geometry.corner(ii) << " (local: " << geometry.local(geometry.corner(ii)) << ")\n";
} // ... printIntersection(...)


/**
 * \brief Checks if intersection contains the given global_point.
 *
 *        Returns true, if global_point lies on the line between the corners of intersection.
 */
template <class G, class I, class D>
  typename std::enable_if< Dune::Intersection< G, I>::dimension == 2, bool >::type
contains(const Dune::Intersection< G, I>& intersection,
         const Dune::FieldVector<D, 2>& global_point,
         const D& tolerance = DSC::FloatCmp::DefaultEpsilon<D>::value())
{
  const auto& geometry = intersection.geometry();
  // get the global coordinates of the intersections corners
  assert(geometry.corners() == 2);
  const auto corner_0 = geometry.corner(0);
  const auto corner_1 = geometry.corner(1);
  // A line is given by $y = a*x + b$. Searching for a and b fails for certain intersections (for instance those
  // parallel to the y axis. So in order to check if the point is on the line between the corners we consider the
  // vectors pointing from the point to each corner. If those are not orthogonal to the intersections normal, the point
  // cannot lie on the line between the two corners.
  const auto normal = intersection.centerUnitOuterNormal();
  for (auto vector : {corner_0 - global_point, corner_1 - global_point})
    if (vector*normal > tolerance)
      return false;
  // Now that we know the point is on the line, check if it is outside the bounding box of the two corners.
  if (   DSC::FloatCmp::lt(global_point[0], std::min(corner_0[0], corner_1[0]), tolerance)
      || DSC::FloatCmp::gt(global_point[0], std::max(corner_0[0], corner_1[0]), tolerance))
    return false;
  if (   DSC::FloatCmp::lt(global_point[1], std::min(corner_0[1], corner_1[1]), tolerance)
      || DSC::FloatCmp::gt(global_point[1], std::max(corner_0[1], corner_1[1]), tolerance))
    return false;
  // At this point we cannot reject the assumtion that the point lies on the line between the two corners.
  return true;
} // ... contains(...)

/**
 * \brief Checks if global_point lies within the plane spanned by the first three corners of intersection.
 * \sa http://math.stackexchange.com/questions/684141/check-if-a-point-is-on-a-plane-minimize-the-use-of-multiplications-and-divisio
 */
template <class G, class I, class D>
  typename std::enable_if< Dune::Intersection< G, I>::dimension == 3, bool >::type
contains(const Dune::Intersection< G, I>& intersection,
         const Dune::FieldVector<D, 3>& global_point,
         const D& tolerance = DSC::FloatCmp::DefaultEpsilon<D>::value())
{
  const auto& geometry = intersection.geometry();
  // get the global coordinates of the intersections corners, there should be at least 3 (ignore the fourth if there is
  // one, 3 points is enough in 3d)
  assert(geometry.corners() >= 3);
  std::vector<Dune::FieldVector<D, 3>> points(4, Dune::FieldVector<D, 3>(0.));
  for (size_t ii = 0; ii < 3; ++ii)
    points[ii] = geometry.corner(ii);
  points[3] = global_point;
  // form a matrix of these points, where the top three entries of each column are given by the three entries of each
  // point and the bottom row is given by one, i.e.
  // a_0 b_0 c_0 d_0
  // a_1 b_1 c_1 d_1
  // a_2 b_2 c_2 d_2
  //   1   1   1   1
  FieldMatrix<D, 4, 4> matrix(1.); // ensures the 1 on the last row
  for (size_t ii = 0; ii < 3; ++ii) // only set the first three rows
    for (size_t jj = 0; jj < 4; ++jj)
      matrix[ii][jj] = points[jj][ii];
  // the point lies on the plane given by the corners if the determinant of this matrix is zero
  const D det = matrix.determinant();
  return std::abs(det) < tolerance;
} // ... contains(...)


/** Check whether a spatial point lies on an intersection.
*
* @param[in] intersection The intersection
* @param[in] globalPoint A Dune::FieldVector with the global coordinates of the point
* @return Returns true if the point lies on the intersection, false otherwise.
*/
template< class IntersectionType, class FieldType, int dim >
bool intersectionContains( const IntersectionType& intersection, const Dune::FieldVector< FieldType, dim >& globalPoint )
{
  // map global coordinates to local coordinates of the intersection
  const auto& intersectionGeometry = intersection.geometry();
  const auto& localPoint = intersectionGeometry.local(globalPoint);

  // get codim 1 reference element
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 3)
  const auto& refElement = ReferenceElements< FieldType, dim-1 >::general(intersectionGeometry.type());
#else
  const auto& refElement = GenericReferenceElements< FieldType, dim-1 >::general(intersectionGeometry.type());
#endif
  // check whether reference element contains the local coordinates
  return refElement.checkInside(localPoint);
} // end function intersectionContains

} // end namespace Grid
} // end of namespace Stuff
} // end namespace Dune

#endif // DUNE_STUFF_GRID_INTERSECTION_HH
