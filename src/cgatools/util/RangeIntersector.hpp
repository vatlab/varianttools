// Copyright 2010 Complete Genomics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License"); you
// may not use this file except in compliance with the License. You
// may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.

#ifndef CGATOOLS_UTIL_RANGE_INTERSECTOR_HPP_
#define CGATOOLS_UTIL_RANGE_INTERSECTOR_HPP_ 1

//! @file RangeIntersector.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"

#include <functional>
#include <map>
#include <vector>
#include <algorithm>

namespace cgatools { namespace util {

    //! Class that maintains a mapping from key to value, where the key
    //! is a type with range semantics, for the purpose of finding the
    //! intersection of a range to a set of ranges.
    //! Requirements:
    //! - TRange class {a, b} with two boundaries a and b of type TBoundary
    //!   with strict total order defined by BoundaryLess operator. Boundaries
    //!   of a given range are obtained using GetBoundary operator that takes
    //!   the range as the first argument and the boundary index (0/1) as the
    //!   second argument.
    //! - Overlap operator such that {a, b} never overlaps {x, y} if b < x.
    //! This implementation uses an RB-tree where every node is augmented with
    //! the maximum upper range boundary for all ranges in the subtree defined
    //! by the node (see CLR book section 15.3, "Interval trees").
    //! The RB-tree itself is a left-leaning RB-tree variant described by
    //! Robert Sedgewick: http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf .
    //! The implementation is adopted from:
    //! http://www.cs.princeton.edu/~rs/talks/LLRB/Java/RedBlackBST.java .
    //! NOTES:
    //! - Currently the tree orders ranges on left boundary and nothing else.
    //!   Since intersect() returns ranges in the tree order, it may be useful
    //!   to allow custom comparators for TRange; however, they'd still need
    //!   to order on the left boundary first for correctness.
    //! - For now I limited the interface to what the old RangeIntersector used
    //!   to provide; feel free to extend it as needed.
    //! - There are no parent pointers, so if we ever need iterators, we will
    //!   need to maintain stacks in them (or add parent pointers, which is
    //!   somewhat costly).
    template < typename TRange,
               typename TBoundary,
               typename TValue,
               typename Overlap,
               typename GetBoundary,
               typename BoundaryLess = std::less<TBoundary> >
    class IntervalTree
    {
    public:
        typedef IntervalTree<TRange, TBoundary, TValue,
                             Overlap, GetBoundary, BoundaryLess> MyType;
        typedef std::pair<const TRange, TValue> value_type;
        typedef const value_type* QueryResultType;

        //! Creates empty tree.
        IntervalTree(const Overlap& overlap = Overlap(),
                     const GetBoundary& getBoundary = GetBoundary(),
                     const BoundaryLess& boundaryLess = BoundaryLess())
            :   overlap_(overlap),
                getBoundary_(getBoundary),
                boundaryLess_(boundaryLess),
                size_(0), root_(0)
        {
        }

        //! Copy constructor. Provides value copy semantics.
        IntervalTree(const MyType& rhs)
            :   overlap_(rhs.overlap_),
                getBoundary_(rhs.getBoundary_),
                boundaryLess_(rhs.boundaryLess_),
                size_(rhs.size_), root_(copyNodes(rhs.root_))
        {
        }

        //! Assignment. Provides value copy semantics.
        const MyType& operator=(const MyType& rhs)
        {
            MyType t(rhs);
            swap(t);
            return *this;
        }

        ~IntervalTree()
        {
            freeNodes(root_);
        }

        //! Add a range to the set.
        void put(const TRange& range, const TValue& value)
        {
            root_ = doInsert(root_, range, value);
        }

        //! Find the intersection of the given range with all the ranges
        //! in the set.
        void intersect(const TRange& range,
                       std::vector< QueryResultType >& result ) const
        {
            result.clear();
            doSearch(root_, range, result);
        }

        //! Number of ranges in the tree.
        size_t size() const { return size_; }

        //! Swaps the internals of this tree with another.
        void swap(MyType& rhs)
        {
            std::swap(overlap_, rhs.overlap_);
            std::swap(getBoundary_, rhs.getBoundary_);
            std::swap(boundaryLess_, rhs.boundaryLess_);
            std::swap(size_, rhs.size_);
            std::swap(root_, rhs.root_);
        }

        //! Clear all data.
        void clear()
        {
            freeNodes(root_);
            root_ = 0;
            size_ = 0;
        }

        //! Debug function: returns tree depth and validates the
        //! tree structure.
        size_t getMaxDepth() const
        {
            size_t depth, sz = 0;
            depth = doGetDepth(root_, sz);
            CGA_ASSERT(size_ == sz);
            CGA_ASSERT(0 == root_ || depth > 0);
            return depth;
        }

    private: // types and constants
        static const bool RED = true;
        static const bool BLACK = false;

        struct Node
        {
            Node(const TRange& k, const TValue& v)
                :   data_(k, v), color_(RED), left_(0), right_(0)
            {
            }

            const TRange& key() const { return data_.first; }

            value_type data_;
            TBoundary rmax_;
            bool color_;
            Node* left_;
            Node* right_;
        };

    private: // data
        Overlap overlap_;
        GetBoundary getBoundary_;
        BoundaryLess boundaryLess_;
        size_t size_;
        Node* root_;

    private: // methods

        bool isRed(const Node* h) const
        {
            if (0 == h)
                return false;
            else
                return h->color_ == RED;
        }

        bool less(const TRange& a, const TRange& b)
        {
            return boundaryLess_( getBoundary_(a, 0), getBoundary_(b, 0) );
        }

        Node* setRangeMax(Node *h)
        {
            h->rmax_ =  getBoundary_(h->key(), 1);
            if (h->left_ && boundaryLess_(h->rmax_, h->left_->rmax_))
                h->rmax_ = h->left_->rmax_;
            if (h->right_ && boundaryLess_(h->rmax_, h->right_->rmax_))
                h->rmax_ = h->right_->rmax_;
            return h;
        }

        void colorFlip(Node* h)
        {
            h->color_ = !h->color_;
            h->left_->color_ = !h->left_->color_;
            h->right_->color_ = !h->right_->color_;
        }

        Node* rotateRight(Node* h)
        {
            Node* x = h->left_;
            h->left_ = x->right_;
            x->right_ = setRangeMax(h);
            x->color_ = x->right_->color_;
            x->right_->color_ = RED;
            return setRangeMax(x);
        }

        Node* rotateLeft(Node* h)
        {
            Node* x = h->right_;
            h->right_ = x->left_;
            x->left_ = setRangeMax(h);
            x->color_ = x->left_->color_;
            x->left_->color_ = RED;
            return setRangeMax(x);
        }

        Node* doInsert(Node* h, const TRange& k, const TValue& v)
        {
            if (!h)
            {
                ++size_;
                return setRangeMax(new Node(k, v));
            }

            if (less(k, h->key()))
                h->left_ = doInsert(h->left_, k, v);
            else
                h->right_ = doInsert(h->right_, k, v);

            if (isRed(h->right_))
                h = rotateLeft(h);

            if (isRed(h->left_) && isRed(h->left_->left_))
                h = rotateRight(h);

            if (isRed(h->left_) && isRed(h->right_))
                colorFlip(h);

            return setRangeMax(h);
        }

        void doSearch(const Node* h, const TRange& k,
                      std::vector< QueryResultType >& result) const
        {
            if (!h)
                return;

            if (h->left_ && !boundaryLess_( h->left_->rmax_, getBoundary_(k, 0) ))
                doSearch(h->left_, k, result);

            if ( boundaryLess_( getBoundary_(k, 1), getBoundary_(h->key(), 0)) )
                return;

            if (overlap_(h->key(), k))
                result.push_back(&h->data_);

            if (h->right_ && !boundaryLess_(h->right_->rmax_, getBoundary_(k, 0)) )
                doSearch(h->right_, k, result);
        }

        size_t doGetDepth(const Node* h, size_t& sz) const
        {
            if (!h)
                return 0;

            sz += 1;

            // Basic sanity
            CGA_ASSERT(!h->right_ || h->left_);
            // Range max. validity
            CGA_ASSERT(!boundaryLess_(h->rmax_, getBoundary_(h->key(), 1)));
            if (h->left_)
                CGA_ASSERT(!boundaryLess_(h->rmax_, h->left_->rmax_));
            if (h->right_)
                CGA_ASSERT(!boundaryLess_(h->rmax_, h->right_->rmax_));
            // 23-tree correctness; note that the color of the root node is bogus
            CGA_ASSERT(!isRed(h->right_));
            CGA_ASSERT(!(isRed(h->left_) && isRed(h->left_->left_)));
            // debugging end

            return 1 + std::max(doGetDepth(h->left_, sz), doGetDepth(h->right_, sz));
        }

        void freeNodes(Node *h)
        {
            if (!h)
                return;
            freeNodes(h->left_);
            freeNodes(h->right_);
            delete h;
        }

        Node* copyNodes(Node *h)
        {
            if (!h)
                return 0;
            Node* t = new Node(*h);
            t->left_ = t->right_ = 0;
            try
            {
                t->left_ = copyNodes(h->left_);
                t->right_ = copyNodes(h->right_);
            }
            catch (const std::exception& )
            {
                freeNodes(t);
                throw;
            }
            return t;
        }
    };

} } // cgatools::util

#endif // CGATOOLS_UTIL_RANGE_INTERSECTOR_HPP_
