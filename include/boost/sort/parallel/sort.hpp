//----------------------------------------------------------------------------
/// @file sort.hpp
/// @brief This file contains the sort functions of the sort library
///
/// @author Copyright (c) 2016 Francisco Jose Tapia (fjtapia@gmail.com )\n
///         Distributed under the Boost Software License, Version 1.0.\n
///         ( See accompanying file LICENSE_1_0.txt or copy at
///           http://www.boost.org/LICENSE_1_0.txt  )
/// @version 0.1
///
/// @remarks
//-----------------------------------------------------------------------------
#ifndef __BOOST_SORT_PARALLEL_SORT_HPP
#define __BOOST_SORT_PARALLEL_SORT_HPP

#include <boost/sort/parallel/detail/select_block_indirect.hpp>
#include <boost/sort/parallel/detail/parallel_stable_sort.hpp>
#include <boost/sort/parallel/detail/util/compare_traits.hpp>
#include <iterator>

namespace boost
{
namespace sort
{
namespace parallel
{
//
//****************************************************************************
//             USING AND DEFINITIONS
//****************************************************************************
using std::iterator_traits;
using detail::less_ptr_no_null;
using detail::util::compare_iter;
using detail::util::enable_if_not_integral;
//
//############################################################################
//                                                                          ##
//                 S O R T , I N D I R E C T _ S O R T                      ##
//                                                                          ##
//############################################################################
//
//-----------------------------------------------------------------------------
//  function : sort
/// @brief non stable sort, based internally in the intro_sort algorithm
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
//-----------------------------------------------------------------------------
template < class Iter_t, typename Compare = compare_iter< Iter_t > >
void sort (Iter_t first, Iter_t last, Compare comp = Compare ( ))
{
    detail::intro_sort (first, last, comp);
};
//
//-----------------------------------------------------------------------------
//  function : indirect_sort
/// @brief : indirect sort. This function sort a vector of iterators to the
///          elements. And after, move the elements, for to be phisically
///          sorted. This is useful with big objects
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
//-----------------------------------------------------------------------------
template < class Iter_t, typename Compare = compare_iter< Iter_t > >
void indirect_sort (Iter_t first, Iter_t last, Compare comp = Compare ( ))
{
    typedef less_ptr_no_null< Iter_t, Compare > compare_ptr;

    std::vector< Iter_t > v_iter;
    detail::create_index (first, last, v_iter);
    detail::intro_sort (v_iter.begin ( ), v_iter.end ( ), compare_ptr (comp));
    detail::sort_index (first, v_iter);
};
//
//############################################################################
//                                                                          ##
//                                                                          ##
//                    P A R A L L E L _ S O R T                             ##
//                                                                          ##
//                                                                          ##
//############################################################################
//
//-----------------------------------------------------------------------------
//  function : parallel_sort
/// @brief this function implement a non stable parallel sort. The number of
///        threads to use is the HW threads of the machine
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
//-----------------------------------------------------------------------------
template < class Iter_t >
void parallel_sort (Iter_t first, Iter_t last)
{
    typedef compare_iter< Iter_t > Compare;
    detail::select_block_indirect (first, last, Compare());
};
//
//-----------------------------------------------------------------------------
//  function : parallel_sort
/// @brief this function implement a non stable parallel sort.
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param nthread : Number of threads to use in the process. When this value
///                  is lower than 2, the sorting is done with 1 thread
//-----------------------------------------------------------------------------
template < class Iter_t >
void parallel_sort (Iter_t first, Iter_t last, uint32_t nthread)
{
    typedef compare_iter< Iter_t > Compare;
    detail::select_block_indirect (first, last, Compare(), nthread);
};
//-----------------------------------------------------------------------------
//  function : parallel_sort
/// @brief this function implement a non stable parallel sort.
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
/// @param nthread : Number of threads to use in the process. When this value
///                  is lower than 2, the sorting is done with 1 thread
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare >
void parallel_sort (Iter_t first, Iter_t last, Compare comp, uint32_t nthread)
{
    detail::select_block_indirect (first, last, comp, nthread);
};
//
//-----------------------------------------------------------------------------
//  function : parallel_sort
/// @brief : non stable parallel sort.
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare,
           enable_if_not_integral< Compare > * = nullptr >
void parallel_sort (Iter_t first, Iter_t last, Compare comp)
{
    detail::select_block_indirect (first, last, comp);
};

//
//############################################################################
//                                                                          ##
//                                                                          ##
//     S T A B L E _ S O R T , I N D I R E C T _ S T A B L E _ S O R T      ##
//                                                                          ##
//                                                                          ##
//############################################################################
//
//-----------------------------------------------------------------------------
//  function : stable_sort
/// @brief this function implement a single thread stable sort
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare = compare_iter< Iter_t > >
void stable_sort (Iter_t first, Iter_t last, Compare comp = Compare ( ))
{
    detail::spin_sort< Iter_t, Compare > (first, last, comp);
};
//
//-----------------------------------------------------------------------------
//  function : indirect_stable_sort
/// @brief : indirect stable sort. This function make a stable sort of a vector
///          of iterators to the elements. And after, move the elements, for
///          to be phisically sorted. This is useful with big objects
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare = compare_iter< Iter_t > >
void indirect_stable_sort (Iter_t first, Iter_t last,
                           Compare comp = Compare ( ))
{
    typedef less_ptr_no_null< Iter_t, Compare > compare_ptr;
    typedef typename std::vector< Iter_t >::iterator iter_ptr;

    std::vector< Iter_t > v_iter;
    detail::create_index (first, last, v_iter);
    detail::spin_sort< iter_ptr, compare_ptr >
      (v_iter.begin ( ), v_iter.end ( ), compare_ptr (comp));
    detail::sort_index (first, v_iter);
};
//
//############################################################################
//                                                                          ##
//                                                                          ##
//            P A R A L L E L _ S T A B L E _ S O R T                       ##
//                                                                          ##
//                                                                          ##
//############################################################################
//
//-----------------------------------------------------------------------------
//  function : parallel_stable_sort
/// @brief : parallel stable sort algorithm.
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
//-----------------------------------------------------------------------------
template < class Iter_t >
void parallel_stable_sort (Iter_t first, Iter_t last)
{
    typedef compare_iter< Iter_t > Compare;
    detail::parallel_stable_sort< Iter_t, Compare > (first, last);
};
//
//-----------------------------------------------------------------------------
//  function : parallel_stable_sort
/// @brief parallel stable sort.
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param nthread : Number of threads to use in the process. When this value
///                  is lower than 2, the sorting is done with 1 thread
//-----------------------------------------------------------------------------
template < class Iter_t >
void parallel_stable_sort (Iter_t first, Iter_t last, uint32_t nthread)
{
    typedef compare_iter< Iter_t > Compare;
    detail::parallel_stable_sort< Iter_t, Compare > (first, last, nthread);
};
//
//-----------------------------------------------------------------------------
//  function : parallel_stable_sort
/// @brief : parallel stable sort.
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare,
           enable_if_not_integral< Compare > * = nullptr >
void parallel_stable_sort (Iter_t first, Iter_t last, Compare comp)
{
    detail::parallel_stable_sort< Iter_t, Compare > (first, last, comp);
};
//
//-----------------------------------------------------------------------------
//  function : parallel_stable_sort
/// @brief parallel stable sort
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
/// @param nthread : Number of threads to use in the process. When this value
///                  is lower than 2, the sorting is done with 1 thread
//-----------------------------------------------------------------------------
template < class Iter_t, typename Compare >
void parallel_stable_sort (Iter_t first, Iter_t last, Compare comp,
                           uint32_t nthread)
{
    detail::parallel_stable_sort< Iter_t, Compare > (first, last, comp,
                                                     nthread);
};
//
//############################################################################
//                                                                          ##
//                                                                          ##
//                       S A M P L E _ S O R T                              ##
//                                                                          ##
//                                                                          ##
//############################################################################
//
//-----------------------------------------------------------------------------
//  function : sample_sort
/// @brief parallel sample sort  algorithm (stable sort)
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
//-----------------------------------------------------------------------------
template < class Iter_t >
void sample_sort (Iter_t first, Iter_t last)
{
    typedef compare_iter< Iter_t > Compare;
    detail::sample_sort< Iter_t, Compare > (first, last);
};
//
//-----------------------------------------------------------------------------
//  function : sample_sort
/// @brief parallel sample sort  algorithm (stable sort)
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param nthread : Number of threads to use in the process. When this value
///                  is lower than 2, the sorting is done with 1 thread
//-----------------------------------------------------------------------------
template < class Iter_t >
void sample_sort (Iter_t first, Iter_t last, uint32_t nthread)
{
    typedef compare_iter< Iter_t > Compare;
    detail::sample_sort< Iter_t, Compare > (first, last, nthread);
};
//
//-----------------------------------------------------------------------------
//  function : sample_sort
/// @brief parallel sample sort  algorithm (stable sort)
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare,
           enable_if_not_integral< Compare > * = nullptr >
void sample_sort (Iter_t first, Iter_t last, Compare comp)
{
    detail::sample_sort< Iter_t, Compare > (first, last, comp);
};
//
//-----------------------------------------------------------------------------
//  function : sample_sort
/// @brief parallel sample sort  algorithm (stable sort)
///
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
/// @param nthread : Number of threads to use in the process. When this value
///                  is lower than 2, the sorting is done with 1 thread
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare >
void sample_sort (Iter_t first, Iter_t last, Compare comp, uint32_t nthread)
{
    detail::sample_sort< Iter_t, Compare > (first, last, comp, nthread);
};
//
//****************************************************************************
};   // End namespace parallel
};   // End namespace sort
};   // End namespace boost
//****************************************************************************
//
#endif
