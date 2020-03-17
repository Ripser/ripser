//----------------------------------------------------------------------------
/// @file low_level.hpp
/// @brief low level functions of create, destroy, move and merge functions
///
/// @author Copyright (c) 2016 Francisco Jose Tapia (fjtapia@gmail.com )\n
///         Distributed under the Boost Software License, Version 1.0.\n
///         ( See accompanying file LICENSE_1_0.txt or copy at
///           http://www.boost.org/LICENSE_1_0.txt  )
/// @version 0.1
///
/// @remarks
//-----------------------------------------------------------------------------
#ifndef __BOOST_SORT_PARALLEL_DETAIL_UTIL_LOW_LEVEL_HPP
#define __BOOST_SORT_PARALLEL_DETAIL_UTIL_LOW_LEVEL_HPP

#include <algorithm>
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
//---------------------------------------------------------------------------
//                         USING SENTENCES
//---------------------------------------------------------------------------
using std::iterator_traits;

//
//-----------------------------------------------------------------------------
//  function : construct
/// @brief create an object in the memory specified by ptr
///
/// @param ptr : pointer to the memory where to create the object
/// @param args : arguments to the constructor
//-----------------------------------------------------------------------------
template < class Value_t, class... Args >
inline void construct (Value_t *ptr, Args &&... args)
{
    (::new (static_cast< void * > (ptr))
         Value_t (std::forward< Args > (args)...));
};
//
//-----------------------------------------------------------------------------
//  function : destroy_object
/// @brief destroy an object in the memory specified by ptr
/// @param ptr : pointer to the object to destroy
//-----------------------------------------------------------------------------
template < class Value_t >
inline void destroy_object (Value_t *ptr)
{
    ptr->~Value_t ( );
};

//
//-----------------------------------------------------------------------------
//  function : init
/// @brief initialize a range of objects with the object val moving across them
///
/// @param first : itertor to the first element to initialize
/// @param last : iterator to the last element to initialize
/// @param val : object used for the initialization
//-----------------------------------------------------------------------------
template < class Iter_t >
void init (Iter_t first, Iter_t last,
           typename iterator_traits< Iter_t >::value_type &val)
{
    if (first == last) return;
    construct (&(*first), std::move (val));

    Iter_t it1 = first, it2 = first + 1;
    while (it2 != last) {
        construct (&(*(it2++)), std::move (*(it1++)));
    };
    val = std::move (*(last - 1));
};
//
//-----------------------------------------------------------------------------
//  function : init_move
/// @brief Move initialized objets
/// @param it_dest : iterator to the final place of the objects
/// @param first : iterator to the first element to move
/// @param last : iterator to the last element to move
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t >
Iter2_t init_move (Iter2_t it_dest, Iter1_t first, Iter1_t last)
{
    while (first != last) *(it_dest++) = std::move (*(first++));
    return it_dest;
};
//
//-----------------------------------------------------------------------------
//  function : uninit_move
/// @brief Move objets to uninitialized memory
///
/// @param ptr : pointer to the memory where to create the objects
/// @param first : iterator to the first element to move
/// @param last : iterator to the last element to move
//-----------------------------------------------------------------------------
template < class Iter_t,
           class Value_t = typename iterator_traits< Iter_t >::value_type >
Value_t *uninit_move (Value_t *ptr, Iter_t first, Iter_t last)
{
    typedef typename iterator_traits< Iter_t >::value_type value2_t;

    static_assert (std::is_same< Value_t, value2_t >::value,
                   "Incompatible iterators\n");

    while (first != last) {
        ::new (static_cast< void * > (ptr++)) Value_t (std::move (*(first++)));
    };
    return ptr;
};
//
//-----------------------------------------------------------------------------
//  function : destroy
/// @brief destroy the elements between first and last
/// @param first : iterator to the first element to destroy
/// @param last : iterator to the last element to destroy
//-----------------------------------------------------------------------------
template < class Iter_t >
void destroy (Iter_t first, const Iter_t last)
{
    typedef typename iterator_traits< Iter_t >::value_type value_t;
    while (first != last) {
        (&(*(first++)))->~value_t ( );
    };
};
//
//-----------------------------------------------------------------------------
//  function : full_merge
/// @brief Merge two contiguous buffers pointed by buf1 and buf2, and put
///        in the buffer pointed by buf_out
///
/// @param buf1 : iterator to the first element in the first buffer
/// @param end_buf1 : final iterator of first buffer
/// @param buf2 : iterator to the first iterator to the second buffer
/// @param end_buf2 : final iterator of the second buffer
/// @param buf_out : buffer where move the elements merged
/// @param comp : comparison object
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
Iter2_t full_merge (Iter1_t buf1, const Iter1_t end_buf1, Iter1_t buf2,
                    const Iter1_t end_buf2, Iter2_t buf_out, Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type value1_t;
    typedef typename iterator_traits< Iter2_t >::value_type value2_t;
    static_assert (std::is_same< value1_t, value2_t >::value,
                   "Incompatible iterators\n");


    while ((buf1 != end_buf1) and (buf2 != end_buf2)) {
        *(buf_out++) = (not comp (*buf2, *buf1)) ? std::move (*(buf1++))
                                                 : std::move (*(buf2++));
    };
    return (buf1 == end_buf1) ? init_move (buf_out, buf2, end_buf2)
                              : init_move (buf_out, buf1, end_buf1);
};
//
//-----------------------------------------------------------------------------
//  function : uninit_full_merge
/// @brief Merge two contiguous buffers pointed by first1 and first2, and put
///        in the uninitialized buffer pointed by it_out
///
/// @param first1 : iterator to the first element in the first buffer
/// @param last1 : last iterator of the first buffer
/// @param first2 : iterator to the first element to the second buffer
/// @param last2 : final iterator of the second buffer
/// @param it_out : uninitialized buffer where move the elements merged
/// @param comp : comparison object
//-----------------------------------------------------------------------------
template < class Iter_t, class Value_t, class Compare >
Value_t *uninit_full_merge (Iter_t first1, const Iter_t last1, Iter_t first2,
                            const Iter_t last2, Value_t *it_out, Compare comp)
{
    typedef typename iterator_traits< Iter_t >::value_type type1;
    static_assert (std::is_same< Value_t, type1 >::value,
                   "Incompatible iterators\n");

    while (first1 != last1 and first2 != last2) {
        construct ((it_out++), (not comp (*first2, *first1))
                                   ? std::move (*(first1++))
                                   : std::move (*(first2++)));
    };
    return (first1 == last1) ? uninit_move (it_out, first2, last2)
                             : uninit_move (it_out, first1, last1);
};
//
//---------------------------------------------------------------------------
//  function : half_merge
/// @brief : Merge two buffers. The first buffer is in a separate memory.
///          The second buffer have a empty space before buf2 of the same size
///          than the (end_buf1 - buf1)
///
/// @param buf1 : iterator to the first element of the first buffer
/// @param end_buf1 : iterator to the last element of the first buffer
/// @param buf2 : iterator to the first element of the second buffer
/// @param end_buf2 : iterator to the last element of the second buffer
/// @param buf_out : iterator to the first element to the buffer where put
///                  the result
/// @param comp : object for Compare two elements of the type pointed
///                by the Iter1_t and Iter2_t
//---------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
Iter2_t half_merge (Iter1_t buf1, const Iter1_t end_buf1, Iter2_t buf2,
                    const Iter2_t end_buf2, Iter2_t buf_out, Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type value1_t;
    typedef typename iterator_traits< Iter2_t >::value_type value2_t;
    static_assert (std::is_same< value1_t, value2_t >::value,
                   "Incompatible iterators\n");

    while ((buf1 != end_buf1) and (buf2 != end_buf2)) {
        *(buf_out++) = (not comp (*buf2, *buf1)) ? std::move (*(buf1++))
                                                 : std::move (*(buf2++));
    };
    return (buf2 == end_buf2) ? init_move (buf_out, buf1, end_buf1) : end_buf2;
};
//
//-----------------------------------------------------------------------------
//  function : in_place_merge_uncontiguous
/// @brief : merge two uncontiguous buffers, placing the results in the buffers
///          Use an auxiliary buffer pointed by aux
///
/// @param src1 : iterator to the first element of the first buffer
/// @param end_src1 : last iterator  of the first buffer
/// @param src2 : iterator to the first element of the second buffer
/// @param end_src2 : last iterator  of the second buffer
/// @param aux  : iterator to the first element of the auxiliary buffer
/// @param comp : object for to Compare elements
/// @return true : not changes done,  false : changes in the buffers
/// @remarks
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Iter3_t, class Compare >
bool in_place_merge_uncontiguous (Iter1_t src1, const Iter1_t end_src1,
                                  Iter2_t src2, const Iter2_t end_src2,
                                  Iter3_t aux, Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    typedef typename iterator_traits< Iter3_t >::value_type type3;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");
    static_assert (std::is_same< type3, type2 >::value,
                   "Incompatible iterators\n");

    if (src1 == end_src1 or src2 == end_src2 or
        not comp (*src2, *(end_src1 - 1)))
        return true;

    while (src1 != end_src1 and not comp (*src2, *src1)) ++src1;

    Iter3_t const end_aux = aux + (end_src1 - src1);
    Iter2_t src2_first = src2;
    init_move (aux, src1, end_src1);

    while ((src1 != end_src1) and (src2 != end_src2)) {
        *(src1++) = std::move ((not comp (*src2, *aux)) ? *(aux++) : *(src2++));
    }

    if (src2 == end_src2) {
        while (src1 != end_src1) *(src1++) = std::move (*(aux++));
        init_move (src2_first, aux, end_aux);
    }
    else
    {
        half_merge (aux, end_aux, src2, end_src2, src2_first, comp);
    };
    return false;
};

//
//-----------------------------------------------------------------------------
//  function : in_place_merge
/// @brief : merge two contiguous buffers,using an auxiliary buffer pointed
///          by buf. The results are in src1 and src2
///
/// @param src1: iterator to the first position of the first buffer
/// @param src2: final iterator of the first buffer and first iterator
///              of the second buffer
/// @param end_src2 : final iterator of the second buffer
/// @param buf  : iterator to buffer used as auxiliary memory
/// @param comp : object for to Compare elements
/// @return true : not changes done,  false : changes in the buffers
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
bool in_place_merge (Iter1_t src1, Iter1_t src2, Iter1_t end_src2, Iter2_t buf,
                     Compare comp)
{
    typedef typename iterator_traits< Iter1_t >::value_type type1;
    typedef typename iterator_traits< Iter2_t >::value_type type2;
    static_assert (std::is_same< type1, type2 >::value,
                   "Incompatible iterators\n");

    if (src1 == src2 or src2 == end_src2 or not comp (*src2, *(src2 - 1)))
        return true;

    Iter1_t end_src1 = src2;
    while (src1 != end_src1 and not comp (*src2, *src1)) ++src1;

    if (src1 == end_src1) return false;

    size_t nx = end_src1 - src1;
    init_move (buf, src1, end_src1);
    half_merge (buf, buf + nx, src2, end_src2, src1, comp);
    return false;
};
//
//****************************************************************************
}; //    End namespace util
}; //    End namespave detail
}; //    End namespace parallel
}; //    End namespace sort
}; //    End namespace boost
//****************************************************************************
//
#endif
