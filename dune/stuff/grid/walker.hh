// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_WALKER_HH
#define DUNE_STUFF_GRID_WALKER_HH

#include <vector>
#include <memory>
#include <type_traits>
#include <functional>

#if 1 // HAVE_TBB
# include <tbb/blocked_range.h>
# include <tbb/parallel_reduce.h>
# include <tbb/tbb_stddef.h>
#endif

#include <dune/stuff/grid/entity.hh>
#include <dune/stuff/grid/intersection.hh>

#include "walker/functors.hh"
#include "walker/apply-on.hh"
#include "walker/wrapper.hh"

namespace Dune {
namespace Stuff {
namespace Grid {


template< class GridViewImp >
class Walker
  : public Functor::Codim0And1< GridViewImp >
{
  typedef Walker< GridViewImp > ThisType;
public:
  typedef GridViewImp GridViewType;
  typedef typename Stuff::Grid::Entity< GridViewType >::Type       EntityType;
  typedef typename Stuff::Grid::Intersection< GridViewType >::Type IntersectionType;

  Walker(const GridViewType& grd_vw)
    : grid_view_(grd_vw)
  {}

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  void add(std::function< void(const EntityType&) > lambda,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    codim0_functors_.emplace_back(new internal::Codim0LambdaWrapper< GridViewType >(lambda, where));
  }

  void add(Functor::Codim0< GridViewType >& functor,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    codim0_functors_.emplace_back(
          new internal::Codim0FunctorWrapper<GridViewType, Functor::Codim0< GridViewType > >(functor, where));
  }

  void add(Functor::Codim1< GridViewType >& functor,
           const ApplyOn::WhichIntersection< GridViewType >* where = new ApplyOn::AllIntersections< GridViewType >())
  {
    codim1_functors_.emplace_back(
          new internal::Codim1FunctorWrapper<GridViewType, Functor::Codim1< GridViewType > >(functor, where));
  }

  void add(Functor::Codim0And1< GridViewType >& functor,
           const ApplyOn::WhichEntity< GridViewType >* which_entities = new ApplyOn::AllEntities< GridViewType >(),
           const ApplyOn::WhichIntersection< GridViewType >* which_intersections
              = new ApplyOn::AllIntersections< GridViewType >())
  {
    codim0_functors_.emplace_back(
          new internal::Codim0FunctorWrapper<GridViewType, Functor::Codim0And1< GridViewType > >(functor,
                                                                                                 which_entities));
    codim1_functors_.emplace_back(
          new internal::Codim1FunctorWrapper<GridViewType, Functor::Codim0And1< GridViewType > >(functor,
                                                                                                 which_intersections));
  }

  void add(Functor::Codim0And1< GridViewType >& functor,
           const ApplyOn::WhichIntersection< GridViewType >* which_intersections,
           const ApplyOn::WhichEntity< GridViewType >* which_entities = new ApplyOn::AllEntities< GridViewType >())
  {
    codim0_functors_.emplace_back(
          new internal::Codim0FunctorWrapper<GridViewType, Functor::Codim0And1< GridViewType > >(functor,
                                                                                                 which_entities));
    codim1_functors_.emplace_back(
          new internal::Codim1FunctorWrapper<GridViewType, Functor::Codim0And1< GridViewType > >(functor,
                                                                                                 which_intersections));
  }

  void add(ThisType& other,
           const ApplyOn::WhichEntity< GridViewType >* which_entities = new ApplyOn::AllEntities< GridViewType >(),
           const ApplyOn::WhichIntersection< GridViewType >* which_intersections
              = new ApplyOn::AllIntersections< GridViewType >())
  {
    if (&other == this)
      DUNE_THROW(Stuff::Exceptions::internal_error, "Do not add a Walker to itself!");
    codim0_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_entities));
    codim1_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_intersections));
  } // ... add(...)

  void add(ThisType& other,
           const ApplyOn::WhichIntersection< GridViewType >* which_intersections,
           const ApplyOn::WhichEntity< GridViewType >* which_entities = new ApplyOn::AllEntities< GridViewType >())
  {
    if (&other == this)
      DUNE_THROW(Stuff::Exceptions::internal_error, "Do not add a Walker to itself!");
    codim0_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_entities));
    codim1_functors_.emplace_back(new internal::WalkerWrapper<GridViewType, ThisType>(other, which_intersections));
  } // ... add(...)

  void clear()
  {
    codim0_functors_.clear();
    codim1_functors_.clear();
  } // ... clear()

  virtual void prepare()
  {
    for (auto& functor : codim0_functors_)
      functor->prepare();
    for (auto& functor : codim1_functors_)
      functor->prepare();
  } // ... prepare()

  bool apply_on(const EntityType& entity) const
  {
    for (const auto& functor : codim0_functors_)
      if (functor->apply_on(grid_view_, entity))
        return true;
    return false;
  } // ... apply_on(...)

  bool apply_on(const IntersectionType& intersection) const
  {
    for (const auto& functor : codim1_functors_)
      if (functor->apply_on(grid_view_, intersection))
        return true;
    return false;
  } // ... apply_on(...)

  virtual void apply_local(const EntityType& entity)
  {
    for (auto& functor : codim0_functors_)
      if (functor->apply_on(grid_view_, entity))
        functor->apply_local(entity);
  } // ... apply_local(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity)
  {
    for (auto& functor : codim1_functors_)
      if (functor->apply_on(grid_view_, intersection))
        functor->apply_local(intersection, inside_entity, outside_entity);
  } // ... apply_local(...)

  virtual void finalize()
  {
    for (auto& functor : codim0_functors_)
      functor->finalize();
    for (auto& functor : codim1_functors_)
      functor->finalize();
  } // ... finalize()

  void walk(const bool clear_stack = true)
  {
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
      walk_range(DSC::viewRange(grid_view_));
    } // only do something, if we have to

    // finalize functors
    finalize();

    // clear the stack of functors
    if (clear_stack) clear();
  } // ... walk(...)

protected:
#if 1 // HAVE_TBB
  template< class PartioningType, class WalkerType >
  struct Body
  {
    Body(WalkerType& walker, PartioningType& partitioning) :
      walker_(walker), partitioning_(partitioning) {}
    Body(Body& other, tbb::split /*split*/)
      : walker_(other.walker_), partitioning_(other.partitioning_) {}

    void operator()(const tbb::blocked_range<std::size_t> &range)
    {
      // for all partitions in tbb-range
      for(std::size_t p = range.begin(); p != range.end(); ++p)
        walker_.walk_range(partitioning_.partition(p));
    }
    void join(Body& /*other*/)
    {}

    WalkerType& walker_;
    PartioningType& partitioning_;
  }; // struct Body

public:
  template< class PartioningType >
  void tbb_walk(PartioningType& partitioning, const bool clear_stack = true)
  {
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
      tbb::blocked_range< std::size_t > range(0, partitioning.partitions());

      Body< PartioningType, ThisType > body(*this, partitioning);
      tbb::parallel_reduce(range, body);
    } // only do something, if we have to

    // finalize functors
    finalize();

    // clear the stack of functors
    if (clear_stack) clear();
  } // ... tbb_walk(...)
#endif

protected:
  template< class EntityRange >
  void walk_range(const EntityRange& entity_range)
  {
    for(const EntityType& entity : entity_range) {
      // apply codim0 functors
      apply_local(entity);

      // only walk the intersections, if there are codim1 functors present
      if (codim1_functors_.size() > 0) {
        // walk the intersections
        const auto intersection_it_end = grid_view_.iend(entity);
        for (auto intersection_it = grid_view_.ibegin(entity);
             intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;

          // apply codim1 functors
          if (intersection.neighbor()) {
            const auto neighbor_ptr = intersection.outside();
            const auto& neighbor = *neighbor_ptr;
            apply_local(intersection, entity, neighbor);
          } else
            apply_local(intersection, entity, entity);

        } // walk the intersections
      } // only walk the intersections, if there are codim1 functors present
    }
  }

  const GridViewType& grid_view_;
  std::vector< std::unique_ptr< internal::Codim0Object<GridViewType> > > codim0_functors_;
  std::vector< std::unique_ptr< internal::Codim1Object<GridViewType> > > codim1_functors_;
}; // class Walker

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_WALKER_HH