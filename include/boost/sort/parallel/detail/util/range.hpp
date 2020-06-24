//----------------------------------------------------------------------------
/// @file range.hpp
/// @brief Define a range [first, last), and the associated operations
///
/// @author Copyright (c) 2016 Francisco José Tapia (fjtapia@gmail.com )\n
///         Distributed under the Boost Software License, Version 1.0.\n
///         ( See accompanyingfile LICENSE_1_0.txt or copy at
///           http://www.boost.org/LICENSE_1_0.txt  )
/// @version 0.1
///
/// @remarks
//-----------------------------------------------------------------------------
#ifndef __BOOST_SORT_PARALLEL_DETAIL_UTIL_RANGE_HPP
#define __BOOST_SORT_PARALLEL_DETAIL_UTIL_RANGE_HPP

#include <boost/sort/parallel/detail/util/low_level.hpp>
#include <cassert>
#include <functional>
#include <iterator>
#include <memory>
#include <type_traits>
#include <vector>

namespace boost
{
namespace sort
{
namespace parallel
{
namespace detail
{
namespace util
{
//----------------------------------------------------------------------------
//               USING SENTENCES
//----------------------------------------------------------------------------
using std::iterator_traits;

///---------------------------------------------------------------------------
/// @struct range
/// @brief this represent a range between two iterators
/// @remarks
//----------------------------------------------------------------------------
template < class Iter_t >
struct range
{
    Iter_t first, last;
    //
    //------------------------------------------------------------------------
    //  function : range
    /// @brief  empty constructor
    //------------------------------------------------------------------------
    range (void){};
    //
    //------------------------------------------------------------------------
    //  function : range
    /// @brief  constructor with two parameters
    /// @param frs : iterator to the first element
    /// @param lst : iterator to the last element
    //-----------------------------------------------------------------------
    range (Iter_t frs, Iter_t lst) : first (frs), last (lst){};
    //
    //-----------------------------------------------------------------------
    //  function : empty
    /// @brief indicate if the range is empty
    /// @return  true : empty false : not empty
    //-----------------------------------------------------------------------
    bool empty (void) const { return (first == last); };
    //
    //-----------------------------------------------------------------------
    //  function : not_empty
    /// @brief indicate if the range is not empty
    /// @return  true : not empty false : empty
    //-----------------------------------------------------------------------
    bool not_empty (void) const { return (first != last); };
    //
    //-----------------------------------------------------------------------
    //  function : valid
    /// @brief  Indicate if the range is well constructed, and valid
    /// @return true : valid,  false : not valid
    //-----------------------------------------------------------------------
    bool valid (void) const { return ((last - first) >= 0); };
    //
    //-----------------------------------------------------------------------
    //  function : size
    /// @brief  return the size of the range
    /// @return size
    //-----------------------------------------------------------------------
    size_t size (void) const { return (last - first); };
    //
    //------------------------------------------------------------------------
    //  function : front
    /// @brief return an iterator to the first element of the range
    /// @return iterator
    //-----------------------------------------------------------------------
    Iter_t front (void) const { return first; };
    //
    //-------------------------------------------------------------------------
    //  function : back
    /// @brief return an iterator to the last element of the range
    /// @return iterator
    //-------------------------------------------------------------------------
    Iter_t back (void) const { return (last - 1); };
};

//
//-----------------------------------------------------------------------------
//  function : concat
/// @brief concatenate two contiguous ranges
/// @param it1 : first range
/// @param it2 : second range
/// @return  range resulting of the concatenation
//-----------------------------------------------------------------------------
template < class Iter_t >
range< Iter_t > concat (const range< Iter_t > &it1, const range< Iter_t > &it2)
{
    return range< Iter_t > (it1.first, it2.last);
};
//
//-----------------------------------------------------------------------------
//  function : init_move
/// @brief Move initialized objets from the range src to dest
/// @param dest : range where move the objects
/// @param src : range from where move the objects
/// @return range with the objects moved and the size adjusted
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t >
inline range< Iter2_t > init_move (const range< Iter2_t > &dest,
                                   const range< Iter1_t > &src)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");

    if (src.size ( ) == 0) {
        return range< Iter2_t > (dest.first, dest.first);
    };
    init_move (dest.first, src.first, src.last);
    return range< Iter2_t > (dest.first, dest.first + src.size ( ));
};
//-----------------------------------------------------------------------------
//  function : uninit_move
/// @brief Move uninitialized objets from the range src creating them in  dest
///
/// @param dest : range where move and create the objects
/// @param src : range from where move the objects
/// @return range with the objects moved and the size adjusted
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t >
inline range< Iter2_t > uninit_move (const range< Iter2_t > &dest,
                                     const range< Iter1_t > &src)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");

    if (src.size ( ) == 0) {
        return range< Iter2_t > (dest.first, dest.first);
    };
    uninit_move (dest.first, src.first, src.last);
    return range< Iter2_t > (dest.first, dest.first + src.size ( ));
};
//
//-----------------------------------------------------------------------------
//  function : destroy
/// @brief destroy a range of objects
/// @param rng : range to destroy
//-----------------------------------------------------------------------------
template < class Iter_t >
inline void destroy (range< Iter_t > rng)
{
    detail::util::destroy (rng.first, rng.last);
};
//
//-----------------------------------------------------------------------------
//  function : init
/// @brief initialize a range of objects with the object val moving across them
/// @param rng : range of elements not initialized
/// @param val : object used for the initialization
/// @return range initialized
//-----------------------------------------------------------------------------
template < class Iter_t >
inline range< Iter_t >
    init (const range< Iter_t > &rng,
          typename iterator_traits< Iter_t >::value_type &val)
{
    init (rng.first, rng.last, val);
    return rng;
};
//
//-----------------------------------------------------------------------------
//  function : is_mergeable
/// @brief : indicate if two ranges have a possible merge
/// @param src1 : first range
/// @param src2 : second range
/// @param comp : object for to compare elements
/// @return true : they can be merged
///         false : they can't be merged
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
inline bool is_mergeable (const range< Iter1_t > &src1,
                          const range< Iter2_t > &src2, Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");

    return comp (*(src2.front ( )), *(src1.back ( )));
};
//
//-----------------------------------------------------------------------------
//  function : full_merge
/// @brief Merge two contiguous ranges src1 and src2, and put the result in
///        the range dest, returning the range merged
///
/// @param dest : range where locate the lements merged. the size of dest
///               must be  greater or equal than the sum of the sizes of
///               src1 and src2
/// @param src1 : first range to merge
/// @param src2 : second range to merge
/// @param comp : comparison object
/// @return range with the elements merged and the size adjusted
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Iter3_t, class Compare >
inline range< Iter3_t > full_merge (const range< Iter3_t > &dest,
                                    const range< Iter1_t > &src1,
                                    const range< Iter2_t > &src2, Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    typedef typename iterator_traits< Iter3_t >::value_type type3;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");
    static_assert (std::is_same< type3, type2 >::value,
                   "Incompatible iterators\n");

    return range< Iter3_t > (dest.first,
                             full_merge (src1.first, src1.last, src2.first,
                                         src2.last, dest.first, comp));
};

//-----------------------------------------------------------------------------
//  function : uninit_full_merge
/// @brief Merge two contiguous uninitialized ranges src1 and src2, and create
///        and move the result in the uninitialized range dest, returning the
///        range merged
//
/// @param dest : range where locate the elements merged. the size of dest
///               must be  greater or equal than the sum of the sizes of
///               src1 and src2. Initially is uninitialize memory
/// @param src1 : first range to merge
/// @param src2 : second range to merge
/// @param comp : comparison object
/// @return range with the elements merged and the size adjusted
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Value_t, class Compare >
inline range< Value_t * > uninit_full_merge (const range< Value_t * > &dest,
                                             const range< Iter1_t > &src1,
                                             const range< Iter2_t > &src2,
                                             Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");
    static_assert (std::is_same< Value_t, type2 >::value,
                   "Incompatible iterators\n");

    return range< Value_t * > (
        dest.first, uninit_full_merge (src1.first, src1.last, src2.first,
                                       src2.last, dest.first, comp));
};
//
//---------------------------------------------------------------------------
//  function : half_merge
/// @brief : Merge two initialized buffers. The first buffer is in a separate
///          memory
//
/// @param dest : range where finish the two buffers merged
/// @param src1 : first range to merge in a separate memory
/// @param src2 : second range to merge, in the final part of the
///               range where deposit the final results
/// @param comp : object for compare two elements of the type pointed
///               by the Iter1_t and Iter2_t
/// @return : range with the two buffers merged
//---------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
inline range< Iter2_t > half_merge (const range< Iter2_t > &dest,
                                    const range< Iter1_t > &src1,
                                    const range< Iter2_t > &src2, Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");

    return range< Iter2_t > (dest.first,
                             half_merge (src1.first, src1.last, src2.first,
                                         src2.last, dest.first, comp));
};
//
//-----------------------------------------------------------------------------
//  function : in_place_merge_uncontiguous
/// @brief : merge two non contiguous ranges src1, src2, using the range
///          aux as auxiliary memory. The results are in the original ranges
//
/// @param src1 : first range to merge
/// @param src2 : second range to merge
/// @param aux : auxiliary range used in the merge
/// @param comp : object for to compare elements
/// @return true : not changes done, false : changes in the buffers
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Iter3_t, class Compare >
inline bool in_place_merge_uncontiguous (const range< Iter1_t > &src1,
                                         const range< Iter2_t > &src2,
                                         const range< Iter3_t > &aux,
                                         Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    typedef typename iterator_traits< Iter3_t >::value_type type3;

    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");
    static_assert (std::is_same< type3, type2 >::value,
                   "Incompatible iterators\n");

    return in_place_merge_uncontiguous (src1.first, src1.last, src2.first,
                                        src2.last, aux.first, comp);
};
//
//-----------------------------------------------------------------------------
//  function : in_place_merge
/// @brief : merge two contiguous ranges ( src1, src2) using buf as
///          auxiliary memory. The results are in the same ranges
/// @param src1 : first range to merge
/// @param src1 : second range to merge
/// @param buf : auxiliary memory used in the merge
/// @param comp : object for to compare elements
/// @return true : not changes done,   false : changes in the buffers
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
inline range< Iter1_t >
    in_place_merge (const range< Iter1_t > &src1, const range< Iter1_t > &src2,
                    const range< Iter2_t > &buf, Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");

    in_place_merge (src1.first, src1.last, src2.last, buf.first, comp);
    return concat (src1, src2);
};

//
//-----------------------------------------------------------------------------
//  function : merge_flow
/// @brief : merge two ranges, as part of a merge the ranges in a list. This
///         function reduce the number of movements compared with inplace_merge
///         when you need to merge a sequence of ranges.
///         This function merge the ranges rbuf and rng2, and the results
///          are in rng1 and rbuf
//
/// @param rng1 : range where locate the first elements of the merge
/// @param rbuf : range which provide the first elements, and where store
///               the last results of the merge
/// @param rng2 : range which provide the last elements to merge
/// @param comp : object for to compare elements
/// @return true : not changes done,  false : changes in the buffers
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
void merge_flow (range< Iter1_t > rng1, range< Iter2_t > rbuf,
                 range< Iter1_t > rng2, Compare cmp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");

    range< Iter2_t > rbx (rbuf);
    range< Iter1_t > rx1 (rng1), rx2 (rng2);
    assert (rbx.size ( ) == rx1.size ( ) and rx1.size ( ) == rx2.size ( ));
    while (rx1.first != rx1.last) {
        *(rx1.first++) = (cmp (*rbx.first, *rx2.first))
                             ? std::move (*(rbx.first++))
                             : std::move (*(rx2.first++));
    };
    if (rx2.first == rx2.last) return;
    if (rbx.first == rbx.last)
        util::init_move (rbuf, rng2);
    else
        util::half_merge (rbuf, rx2, rbx, cmp);
};

//****************************************************************************
}; //    End namespace util
}; //    End namespace detail
}; //    End namespace parallel
}; //    End namespace sort
}; //    End namespace boost
//****************************************************************************
//
#endif
