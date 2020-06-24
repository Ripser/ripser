//----------------------------------------------------------------------------
/// @file intro_sort.hpp
/// @brief Intro Sort algorithm
///
/// @author Copyright (c) 2016 Francisco José Tapia (fjtapia@gmail.com )\n
///         Distributed under the Boost Software License, Version 1.0.\n
///         ( See accompanying file LICENSE_1_0.txt or copy at
///           http://www.boost.org/LICENSE_1_0.txt  )
/// @version 0.1
///
/// @remarks
//-----------------------------------------------------------------------------
#ifndef __BOOST_SORT_PARALLEL_ALGORITHM_INTRO_SORT_HPP
#define __BOOST_SORT_PARALLEL_ALGORITHM_INTRO_SORT_HPP

#include <algorithm>
#include <boost/sort/parallel/detail/heap_sort.hpp>
#include <boost/sort/parallel/detail/indirect.hpp>
#include <boost/sort/parallel/detail/insertion_sort.hpp>
#include <boost/sort/parallel/detail/util/compare_traits.hpp>
#include <boost/sort/parallel/detail/util/nbits.hpp>
#include <iterator>
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
//---------------------------------------------------------------------------
//                   USING SENTENCES
//---------------------------------------------------------------------------
using util::compare_iter;
using util::nbits64;
//
//-----------------------------------------------------------------------------
//  function : mid3
/// @brief : return the iterator to the mid value of the three values passsed
///          as parameters
//
/// @param iter_1 : iterator to the first value
/// @param iter_2 : iterator to the second value
/// @param iter_3 : iterator to the third value
/// @param comp : object for to compare two values
/// @return iterator to mid value
//-----------------------------------------------------------------------------
template < typename Iter_t, typename Compare >
inline Iter_t mid3 (Iter_t iter_1, Iter_t iter_2, Iter_t iter_3, Compare comp)
{
    return comp (*iter_1, *iter_2)
               ? (comp (*iter_2, *iter_3)
                      ? iter_2
                      : (comp (*iter_1, *iter_3) ? iter_3 : iter_1))
               : (comp (*iter_3, *iter_2)
                      ? iter_2
                      : (comp (*iter_3, *iter_1) ? iter_3 : iter_1));
};
//
//-----------------------------------------------------------------------------
//  function : pivot3
/// @brief : receive a range between first and last, calcule the mid iterator
///          with the first, the previous to the last, and the central
///          position. With this mid iterator swap with the first position
//
/// @param first : iterator to the first element
/// @param last : iterator to the last element
/// @param comp : object for to compare two elements
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare >
inline void pivot3 (Iter_t first, Iter_t last, Compare comp)
{
    auto N2 = (last - first) >> 1;
    Iter_t it_val = mid3 (first + 1, first + N2, last - 1, comp);
    std::swap (*first, *it_val);
};
//
//-----------------------------------------------------------------------------
//  function : mid9
/// @brief : return the iterator to the mid value of the nine values passsed
///          as parameters
//
/// @param iter_1 : iterator to the first value
/// @param iter_2 : iterator to the second value
/// @param iter_3 : iterator to the third value
/// @param iter_4 : iterator to the fourth value
/// @param iter_5 : iterator to the fifth value
/// @param iter_6 : iterator to the sixth value
/// @param iter_7 : iterator to the seventh value
/// @param iter_8 : iterator to the eighth value
/// @param iter_9 : iterator to the ninth value
/// @return iterator to the mid value
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare >
inline Iter_t mid9 (Iter_t iter_1, Iter_t iter_2, Iter_t iter_3, Iter_t iter_4,
                    Iter_t iter_5, Iter_t iter_6, Iter_t iter_7, Iter_t iter_8,
                    Iter_t iter_9, Compare comp)
{
    return mid3 (mid3 (iter_1, iter_2, iter_3, comp),
                 mid3 (iter_4, iter_5, iter_6, comp),
                 mid3 (iter_7, iter_8, iter_9, comp), comp);
};
//
//-----------------------------------------------------------------------------
//  function : pivot9
/// @brief : receive a range between first and last, obtain 9 values between
///          the elements  including the first and the previous to the last.
///          Obtain the iterator to the mid value and swap with the first
///          position
//
/// @param first : iterator to the first element
/// @param last : iterator to the last element
/// @param comp : object for to compare two elements
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare >
inline void pivot9 (Iter_t first, Iter_t last, Compare comp)
{
    size_t cupo = (last - first) >> 3;
    Iter_t itaux = mid9 (first + 1, first + cupo, first + 2 * cupo,
                         first + 3 * cupo, first + 4 * cupo, first + 5 * cupo,
                         first + 6 * cupo, first + 7 * cupo, last - 1, comp);
    std::swap (*first, *itaux);
};
//
//-----------------------------------------------------------------------------
//  function : intro_sort_internal
/// @brief : internal function for to divide and sort the ranges
//
/// @param first : iterator to the first element
/// @param last : iterator to the element after the last in the range
/// @param level : level of depth from the top level call
/// @param comp : object for to Compare elements
//-----------------------------------------------------------------------------
template < class Iter_t, typename Compare >
void intro_sort_internal (Iter_t first, Iter_t last, uint32_t level,
                          Compare comp)
{
    typedef typename std::iterator_traits< Iter_t >::value_type value_t;

    const uint32_t nmin = 32;
    size_t nelem = last - first;
    if (nelem < nmin) return insertion_sort (first, last, comp);

    if (level == 0) {
        heap_sort< Iter_t, Compare > (first, last, comp);
        return;
    };

    pivot3 (first, last, comp);

    const value_t &val = const_cast< value_t & > (*first);
    Iter_t c_first = first + 1, c_last = last - 1;

    while (comp (*c_first, val)) ++c_first;
    while (comp (val, *c_last)) --c_last;

    while (not(c_first > c_last)) {
        std::swap (*(c_first++), *(c_last--));
        while (comp (*c_first, val)) ++c_first;
        while (comp (val, *c_last)) --c_last;
    }; // End while

    std::swap (*first, *c_last);
    intro_sort_internal (first, c_last, level - 1, comp);
    intro_sort_internal (c_first, last, level - 1, comp);
};

//
//-----------------------------------------------------------------------------
//  function : intro_sort
/// @brief : function for to sort range [first, last )
/// @param first : iterator to the first element
/// @param last : iterator to the element after the last in the range
/// @param comp : object for to compare elements
//-----------------------------------------------------------------------------
template < class Iter_t, typename Compare >
void intro_sort (Iter_t first, Iter_t last, Compare comp)
{
    auto nelem = last - first;
    assert (nelem >= 0);

    //------------------- check if sort --------------------------------------
    bool sw = true;
    for (Iter_t it1 = first, it2 = first + 1;
         it2 != last and (sw = not comp (*it2, *it1)); it1 = it2++)
        ;
    if (sw) return;

    uint32_t level = ((nbits64 (nelem) - 4) * 3) / 2;
    intro_sort_internal (first, last, level, comp);
};
//
//****************************************************************************
}; //    End namespace detail
}; //    End namespace parallel
}; //    End namespace sort
}; //    End namespace boost
//****************************************************************************
//
#endif
